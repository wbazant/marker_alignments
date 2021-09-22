#!/bin/bash
# Reports multiple mamches split by taxon pairs
# Prints pairs of taxa that match in the sate query, with better / equal / worse identity
if [ ! "$1" ] ; then
  echo "Usage: $0 database.sql"
  exit 1
fi

sqlite3 -header -separator $'\t' "$1" <<< "
select at, bt, sum(is_am_higher) as num_markers_at_higher, sum(is_bm_higher) as num_markers_bm_higher
  from
  (
 select 
       at,
       bt,
       am,
       bm,
       sum(am_higher_identity_than_bm) > sum(bm_higher_identity_than_am) as is_am_higher,
       sum(am_higher_identity_than_bm) < sum(bm_higher_identity_than_am) as is_bm_higher
from   (select 
               a.taxon at,
               b.taxon bt,
               a.marker am,
               b.marker bm,
               max(a.identity) - max(b.identity) between -1e-6 and 1e-6 as am_equal_identity_as_bm,
               max(a.identity) - max(b.identity) > 1e-6 as am_higher_identity_than_bm,
               max(b.identity) - max(a.identity) > 1e-6 as bm_higher_identity_than_am
        from   alignment a,
               alignment b
        where  a.query = b.query
        and a.taxon != b.taxon
        and a.marker != b.marker
        group by a.query, at, bt, am, bm
               )
group  by at, bt, am, bm
) group by at, bt;
"
