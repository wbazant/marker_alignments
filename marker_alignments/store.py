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

marker_query_template='''
  select taxon, marker, {} from (
    select
      a.query,
      a.taxon,
      a.marker,
      {}
    from
      alignment a join (
      select query, sum(weight) as total_weight_for_query
        from alignment group by query
      ) as m
    where a.query = m.query
    group by a.taxon, a.marker, a.query
  ) group by taxon, marker
'''
s_cov="sum(coverage) as marker_coverage"
p_cov="a.coverage * a.weight / (m.total_weight_for_query) as coverage"
s_mrc="sum(weight_fraction) as marker_read_count"
p_wf="a.weight / (m.total_weight_for_query) as weight_fraction"
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
a_tnm = "count(marker) as taxon_num_markers, sum(marker_read_count) as taxon_num_reads"
a_cpm = "avg(marker_coverage) / (?) * 1000000 as cpm"
taxon_coverage_query=taxon_query_template.format(a_cov, marker_coverage_query)
taxon_read_and_marker_count_query=taxon_query_template.format(a_tnm, marker_read_count_query)
taxon_cpm_query=taxon_query_template.format(a_cpm, marker_coverage_query)

taxon_joined_query_template='''
  select taxon, {}
  from ({}) join taxon_num_markers_ref using(taxon)
  group by taxon
'''
a_mar = "num_markers_ref as taxon_total_markers, (sum(marker_read_count * marker_read_count) / num_markers_ref - sum(marker_read_count) * sum(marker_read_count) / (num_markers_ref * num_markers_ref)) as taxon_variance"
taxon_markers_query = taxon_joined_query_template.format(a_mar, marker_all_query)
taxon_all_query=taxon_joined_query_template.format(", ".join([a_cov, a_cpm, a_tnm, a_mar]), marker_all_query)

class AlignmentStore(SqliteStore):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.connect()
        self.do('''create table alignment (
              taxon text not null,
              marker text not null,
              query text not null,
              weight real not null,
              coverage real not null
            );''')
        self.do('''create table taxon_num_markers_ref (
              taxon text not null,
              num_markers_ref integer not null
            );''')
        
    def add_alignment(self, taxon, marker, query, weight, coverage):
        self.do('insert into alignment (taxon, marker, query, weight, coverage) values (?,?,?,?,?)', [ taxon, marker, query, weight, coverage])

    def add_taxon_num_markers_ref(self, taxon, num_markers):
        self.do('insert into taxon_num_markers_ref (taxon, num_markers_ref) values (?,?)', [taxon, num_markers])

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

    def as_taxon_markers(self):
        return self.query(taxon_markers_query)

    def as_taxon_cpm(self,total_reads):
        return self.query(taxon_cpm_query, [total_reads])

    def as_taxon_all(self,total_reads):
        return self.query(taxon_all_query, [total_reads, total_reads])
