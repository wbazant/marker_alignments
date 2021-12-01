

# populated when the module loads
num_reads_arg_count_in_sql = {}
sqls = {}

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

sqls['marker_coverage'] = marker_query_template.format(s_cov, p_cov)
sqls['marker_read_count'] = marker_query_template.format(s_mrc, p_wf)
sqls['marker_cpm'] = marker_query_template.format(s_cpm, p_cov)
num_reads_arg_count_in_sql['marker_cpm'] = 1
sqls['marker_all'] = marker_query_template.format(", ".join([s_cov, s_cpm, s_mrc]), ", ".join([p_cov, p_wf]))
num_reads_arg_count_in_sql['marker_all'] = 1



taxon_query_template='''
  select taxon, {}
  from ({})
  group by taxon
'''

a_cov = "avg(marker_coverage) as coverage"
a_tnm = "sum(marker_read_count) as taxon_num_reads, count(marker) as taxon_num_markers, max(marker_read_count) as taxon_max_reads_in_marker"
a_cpm = "avg(marker_coverage) / (?) * 1000000 as cpm"
sqls['taxon_coverage'] = taxon_query_template.format(a_cov, sqls['marker_coverage'])
sqls['taxon_read_and_marker_count'] = taxon_query_template.format(a_tnm, sqls['marker_read_count'])
sqls['taxon_cpm'] = taxon_query_template.format(a_cpm, sqls['marker_coverage'])
num_reads_arg_count_in_sql['taxon_cpm'] = 1
sqls['taxon_all'] = taxon_query_template.format(", ".join([a_cov, a_cpm, a_tnm]), sqls['marker_all'])
num_reads_arg_count_in_sql['taxon_all'] = 2


output_type_options = [k for k in sqls]

def get_output(alignment_store, output_type, num_reads):
    args = [sqls[output_type]]
    if output_type in num_reads_arg_count_in_sql:
        args.append([num_reads] * num_reads_arg_count_in_sql[output_type])
    return alignment_store.report(*args)

field_formats = {
  "taxon" : "",
  "marker": "",
  "marker_cpm": ":.6f",
  "marker_coverage": ":.6f",
  "marker_read_count": ":.2f",
  "marker_avg_identity": ":.6f",
  "cpm" : ":.6f",
  "coverage" : ":.6f",
  "taxon_num_reads": ":.6f",
  "taxon_num_markers": ":d",
  "taxon_max_reads_in_marker": ":.6f",
}

def write(alignment_store, output_type, output_path, num_reads):
    header, lines = get_output(alignment_store, output_type, num_reads)
    formatter="\t".join(['{' + field_formats[field] +'}' for field in header]) + "\n"
    with open(output_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for line in lines:
            f.write(formatter.format(*line))
