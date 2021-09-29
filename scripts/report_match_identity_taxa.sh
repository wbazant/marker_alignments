#! /bin/bash
if [ ! "$1" ] ; then
  echo "Usage: $0 database.sql"
  exit 1
fi
sqlite3 -header -separator $'\t' "$1" <<< "
select taxon,
  avg(top_identity)
from   (select a.taxon,
        max(a.identity) as top_identity
        from alignment a
        group by a.query, a.taxon
        )
group  by taxon;
"

