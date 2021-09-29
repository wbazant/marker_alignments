#! /bin/bash
if [ ! "$1" ] ; then
  echo "Usage: $0 database.sql"
  exit 1
fi
sqlite3 -header -separator $'\t' "$1" <<< "
select taxon,
  count(*) as num_markers,
  avg(avg_identity) as avg_marker_identity,
  sum(higher_identity) as num_markers_at_least_cluster_average,
  sum(lower_identity) as num_markers_below_cluster_average
from (
  select t1.*, t2.avg_cluster_identity, t2.num_taxa, t2.avg_cluster_identity <= avg_identity as higher_identity, t2.avg_cluster_identity > avg_identity as lower_identity
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
order by -num_markers
"
