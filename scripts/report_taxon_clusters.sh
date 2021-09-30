

#! /bin/bash
if [ ! "$1" ] ; then
  echo "Usage: $0 database.sql"
  exit 1
fi
sqlite3 -header -separator $'\t' "$1" <<< "
select tc.id, tc.taxon, count(distinct a.marker) as num_markers, count(distinct a.query) as num_reads, avg(a.identity) as avg_identity
from taxon_cluster tc, alignment a
where tc.taxon = a.taxon
group by tc.id, tc.taxon
"
