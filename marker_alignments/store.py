import sqlite3

class SqliteStore:
    def __init__(self, db_path = None):
        self.__db_path = db_path
        self.__conn = None
        self.__is_within_transaction = False
        self.__stateful_ops_in_bulk_write = None

    def connect(self):
        if self.__conn:
            return
        self.__conn = sqlite3.connect(self.__db_path or ":memory:", isolation_level=None)

    def do(self, *args):
        self.__conn.execute(*args)
        if self.__is_within_transaction:
            self.__stateful_ops_in_bulk_write +=1
            if self.__stateful_ops_in_bulk_write % 100000 == 0:
                self.__conn.execute("commit transaction")
                self.__conn.execute("begin transaction")

    def query(self, *args):
        return self.__conn.execute(*args)

    def start_bulk_write(self):
        self.__is_within_transaction = True
        self.__stateful_ops_in_bulk_write = 0
        self.__conn.execute("begin transaction")

    def end_bulk_write(self):
        self.__conn.execute("commit transaction")
        self.__is_within_transaction = False
        self.__stateful_ops_in_bulk_write = None


# when splitting read stats by query, do it proportionally to the second power of match identity
# if there are multiple matches in a query + taxon + marker, return identity as max and coverage as weighted average
marker_query_template='''
  select taxon, marker, {} from (
    select
      a.query,
      a.taxon,
      a.marker,
      max(a.identity) as identity,
      {}
    from
      alignment a join (
      select query, sum(identity * identity) as total_weight_for_query
        from alignment group by query
      ) as m
    where a.query = m.query
    group by a.taxon, a.marker, a.query
  ) group by taxon, marker
'''
s_cov="sum(coverage) as marker_coverage"
p_cov="sum(a.coverage * a.identity * a.identity) / (m.total_weight_for_query) as coverage"
s_mrc="sum(weight_fraction) as marker_read_count, avg(identity) as marker_avg_identity"
p_wf="sum(a.identity * a.identity) / (m.total_weight_for_query) as weight_fraction"
s_cpm="sum(coverage) / (?) * 1000000 as marker_cpm"

marker_coverage_query=marker_query_template.format(s_cov, p_cov)
marker_read_count_query=marker_query_template.format(s_mrc, p_wf)
marker_cpm_query=marker_query_template.format(s_cpm, p_cov)
marker_all_query=marker_query_template.format(", ".join([s_cov, s_cpm, s_mrc]), ", ".join([p_cov, p_wf]))

taxon_query_template='''
  select taxon, {}
  from ({})
  group by taxon
'''

a_cov = "avg(marker_coverage) as coverage"
a_tnm = "sum(marker_read_count) as taxon_num_reads, count(marker) as taxon_num_markers, max(marker_read_count) as taxon_max_reads_in_marker"
a_cpm = "avg(marker_coverage) / (?) * 1000000 as cpm"
taxon_coverage_query=taxon_query_template.format(a_cov, marker_coverage_query)
taxon_read_and_marker_count_query=taxon_query_template.format(a_tnm, marker_read_count_query)
taxon_cpm_query=taxon_query_template.format(a_cpm, marker_coverage_query)
taxon_all_query=taxon_query_template.format(", ".join([a_cov, a_cpm, a_tnm]), marker_all_query)

filter_taxa_on_multiple_matches_query = '''
  select a.* from alignment a,
  (
    select taxon,
       count(*) as num_matches,
       sum(is_unique) as num_unique_matches,
       sum(is_best) as num_best_matches,
       sum(is_inferior) as num_inferior_matches
    from   (select a.taxon,
      s.num_taxa == 1 as is_unique,
      s.num_taxa > 1 and s.top_identity - max(a.identity) < 1e-6 as is_best,
      s.num_taxa > 1 and s.top_identity - max(a.identity) > 1e-6 as is_inferior
      from   alignment a,
           (select query,
               Max(identity) as top_identity,
               count(distinct taxon) as num_taxa
          from   alignment
          group  by query) s
      where  a.query = s.query
      group by a.query, a.taxon
      )
    group  by taxon
  ) t
  where a.taxon = t.taxon and (t.num_unique_matches + t.num_best_matches ) > (?) * t.num_matches
'''

class AlignmentStore(SqliteStore):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.connect()
        self.do('''create table alignment (
              taxon text not null,
              marker text not null,
              query text not null,
              identity real not null,
              coverage real not null
            );''')
        
    def add_alignment(self, taxon, marker, query, identity, coverage):
        self.do('insert into alignment (taxon, marker, query, identity, coverage) values (?,?,?,?,?)', [ taxon, marker, query, identity, coverage])

    def modify_table_filter_taxa_on_multiple_matches(self, min_fraction_primary_matches):
        self.query("begin transaction")
        self.query("create table new as " + filter_taxa_on_multiple_matches_query, [min_fraction_primary_matches])
        self.query("alter table alignment rename to alignment_pre_filter_taxa_on_multiple_matches")
        self.query("alter table new rename to alignment")
        self.query("commit transaction")


    def as_marker_coverage(self):
        return self.query(marker_coverage_query)

    def as_marker_read_count(self):
        return self.query(marker_read_count_query)

    def as_marker_cpm(self, total_reads):
        return self.query(marker_cpm_query, [total_reads])

    def as_marker_all(self, total_reads):
        return self.query(marker_all_query, [total_reads])

    def as_taxon_coverage(self):
        return self.query(taxon_coverage_query)

    def as_taxon_read_and_marker_count(self):
        return self.query(taxon_read_and_marker_count_query)

    def as_taxon_cpm(self,total_reads):
        return self.query(taxon_cpm_query, [total_reads])

    def as_taxon_all(self,total_reads):
        return self.query(taxon_all_query, [total_reads, total_reads])
