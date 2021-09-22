#! /bin/bash
# Builds a report from a database based on queries with multiple matches
# Summarizes taxon + num {unique, best, inferior} matches
# Matches are compared on % identity score
if [ ! "$1" ] ; then
  echo "Usage: $0 database.sql"
  exit 1
fi
sqlite3 -header -separator $'\t' "$1" <<< "
select taxon,
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
group  by taxon;
"

