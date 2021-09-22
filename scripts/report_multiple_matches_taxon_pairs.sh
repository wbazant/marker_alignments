#!/bin/bash
# Reports multiple matches split by taxon pairs
# Prints pairs of taxa that match in the same query, with better / equal / worse identity
if [ ! "$1" ] ; then
  echo "Usage: $0 database.sql"
  exit 1
fi

sqlite3 -header -separator $'\t' "$1" <<< "
 select at,
       bt,
       sum(at_higher_identity_than_bt) as num_at_higher,
       sum(at_equal_identity_as_bt) as num_equal_identity,
       sum(bt_higher_identity_than_at) as num_bt_higher
from   (select a.taxon at,
               b.taxon bt,
               max(a.identity) - max(b.identity) between -1e-6 and 1e-6 as at_equal_identity_as_bt,
               max(a.identity) - max(b.identity) > 1e-6 as at_higher_identity_than_bt,
               max(b.identity) - max(a.identity) > 1e-6 as bt_higher_identity_than_at
        from   alignment a,
               alignment b
        where  a.query = b.query
        and a.taxon != b.taxon
        group by a.query, at, bt
               )
group  by at,
          bt
order  by at;
"
