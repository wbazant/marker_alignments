import sqlite3

from marker_alignments.mcl import clusters

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
  where a.taxon = t.taxon and (t.num_unique_matches + t.num_best_matches ) >= (?) * t.num_matches
'''

# markers with only inferior alignments don't count
filter_taxa_on_num_markers_and_reads_query = '''
  select a.* from alignment a,
  (
    select taxon,
    count(distinct marker) as num_markers,
    count(distinct query) as num_reads
    from   (
      select a.taxon, a.marker, a.query
      from   alignment a,
           (select query,
               Max(identity) as top_identity,
               count(distinct taxon) as num_taxa
          from   alignment
          group  by query) s
      where  a.query = s.query
      group by a.query, a.taxon, a.marker
      having s.top_identity - max(a.identity) < 1e-6
      )
    group  by taxon
  ) t
  where a.taxon = t.taxon and t.num_markers >= (?) and t.num_reads >= (?)
'''

filter_taxa_on_avg_identity_query = '''
  select a.* from alignment a,
  (
    select taxon,
      avg(top_identity) as avg_identity
    from   (select a.taxon,
            max(a.identity) as top_identity
            from alignment a
            group by a.query, a.taxon
            )
    group  by taxon
  ) t
  where a.taxon = t.taxon and t.avg_identity >= (?)
'''

filter_taxa_on_cluster_averages_query = '''
  select a.* from alignment a,
  (
    select taxon,
      sum(higher_identity) as num_markers_at_least_cluster_average,
      sum(lower_identity) as num_markers_below_cluster_average
    from (
  select t1.*,
        t2.avg_cluster_identity,
        t2.num_taxa,
        t2.avg_cluster_identity - avg_identity < 1e-6 as higher_identity,
        t2.avg_cluster_identity - avg_identity >= 1e-6 as lower_identity
        from (
          select id, mc.taxon, mc.marker, count(distinct query) as num_matches, avg(identity) as avg_identity
            from marker_cluster mc, alignment a
            where mc.taxon = a.taxon and mc.marker = a.marker
            group by id, mc.taxon, mc.marker
        ) t1, (
        select id, avg(identity) as avg_cluster_identity, count(distinct mc.taxon) as num_taxa
            from marker_cluster mc, alignment a
            where mc.taxon = a.taxon and mc.marker = a.marker
            group by id
        ) t2
        where t1.id = t2.id
       )
    group by taxon
  ) t
  where a.taxon = t.taxon and num_markers_at_least_cluster_average >= (?) * num_markers_below_cluster_average
'''

counts_of_common_matches_in_markers_query = '''
select 
       a.taxon at,
       a.marker am,
       b.taxon bt,
       b.marker bm,
       count(distinct a.query)
from   alignment a,
       alignment b
where  a.query = b.query
group by at, bt, am, bm;
'''

counts_of_common_matches_in_taxa_query = '''
select aa.at, aa.bt, cast (sum_shared  as real) / aaa.num_queries from
(
    select at, bt, count(*) as sum_shared
    from (
      select
           a.taxon at,
           b.taxon bt,
           a.query
      from   alignment a,
           alignment b
      where  a.query = b.query
      group by at, bt, a.query
    ) group by at, bt
) aa,
(
  select taxon, count(distinct query) as num_queries from alignment
  group by taxon
) aaa
where aa.at = aaa.taxon
'''


transform_taxa_on_thresholds_and_clusters_query = '''
select t.mapped_taxon as taxon, a.marker, a.query, a.identity, a.coverage from alignment a,
(
    select tc.taxon as original_taxon, tc.taxon as mapped_taxon
    from taxon_cluster tc, alignment al
    where tc.taxon = al.taxon 
    group by tc.id, tc.taxon
    having avg(al.identity) >= (?)

    union

    select 
      tc.taxon as original_taxon,
      m.mapped_taxon
    from taxon_cluster tc,
    (
      select id, '?' || group_concat(taxon) as mapped_taxon
      from (
        select tc.id,
          tc.taxon,
          count(distinct al.marker) as num_markers,
          count(distinct al.query) as num_reads,
          avg(al.identity) as avg_identity
        from taxon_cluster tc, alignment al
        where tc.taxon = al.taxon 
        group by tc.id, tc.taxon
        having avg(al.identity) < (?)
      ) group by id
      having
      (?) > 0 and count(distinct taxon) >= (?) and sum(num_markers) >= (?) and sum(num_reads) >= (?)
    ) m
    where tc.id = m.id
) t
where a.taxon = t.original_taxon
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


    def _modify_table(self, op, select_query, *args):
        self.query("begin transaction")
        self.query("create table new as " + select_query, *args)
        self.query("alter table alignment rename to alignment_pre_filter_on_" + op)
        self.query("alter table new rename to alignment")
        self.query("commit transaction")

    def modify_table_filter_taxa_on_multiple_matches(self, min_fraction_primary_matches):
        self._modify_table('multiple_matches', filter_taxa_on_multiple_matches_query, [min_fraction_primary_matches])

    def modify_table_filter_taxa_on_num_markers_and_reads(self, min_num_markers, min_num_reads):
        self._modify_table('num_markers', filter_taxa_on_num_markers_and_reads_query, [min_num_markers, min_num_reads])

    def modify_table_filter_taxa_on_avg_identity(self, min_avg_identity):
        self._modify_table('avg_identity', filter_taxa_on_avg_identity_query, [min_avg_identity])

    def modify_table_filter_taxa_on_cluster_averages(self, min_better_cluster_averages_ratio):
        self._modify_table('cluster_averages', filter_taxa_on_cluster_averages_query, [min_better_cluster_averages_ratio])

    def modify_table_transform_taxa_on_thresholds_and_clusters(self, threshold_identity, min_num_taxa_below_identity, min_num_markers_below_identity, min_num_reads_below_identity):
        try_return_unknown_taxa = 1 if min_num_taxa_below_identity or min_num_markers_below_identity or min_num_reads_below_identity else 0
        self._modify_table('thresholds_and_clusters', transform_taxa_on_thresholds_and_clusters_query,
	  [threshold_identity, threshold_identity, try_return_unknown_taxa, min_num_taxa_below_identity, min_num_markers_below_identity, min_num_reads_below_identity])

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

    def cluster_markers_by_matches(self):
        triples = [(at + "\t" + am, bt + "\t" + bm, v)  for at, am, bt, bm, v in self.query(counts_of_common_matches_in_markers_query)]

        result = clusters(triples)

        self.do('''
            create table marker_cluster (
              id number not null,
              taxon text not null,
              marker text not null
            );''')

        self.start_bulk_write()
        for ix in range(0, len(result)):
            cluster = result[ix]
            cluster_id = ix + 1
            for x in cluster:
                taxon, marker = x.split("\t")
                self.do('insert into marker_cluster (id, taxon, marker) values (?,?,?)', [ cluster_id, taxon, marker])
        self.end_bulk_write()

    def cluster_taxa_by_matches(self):
        triples = [(at, bt, v)  for at, bt, v in self.query(counts_of_common_matches_in_taxa_query)]

        result = clusters(triples)

        self.do('''
            create table taxon_cluster (
              id number not null,
              taxon text not null
            );''')

        self.start_bulk_write()
        for ix in range(0, len(result)):
            cluster = result[ix]
            cluster_id = ix + 1
            for taxon in cluster:
                self.do('insert into taxon_cluster (id, taxon) values (?,?)', [ cluster_id, taxon])
        self.end_bulk_write()
