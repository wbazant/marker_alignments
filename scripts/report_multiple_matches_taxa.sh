#! /bin/bash
# Builds a report from a database based on queries with multiple matches
# Summarizes taxon + num {unique, best, inferior} matches
# Matches are compared on % identity score
if [ ! "$1" ] ; then
  echo "Usage: $0 database.sql"
  exit 1
fi
sqlite3 -separator $'\t' "$1" <<< "
select taxon,
       match_type,
       count(*)
from   (select a.taxon,
               CASE
                 WHEN s.num_taxa == 1 THEN 'unique'
                 WHEN s.top_identity - max(a.identity) < 1e-6 THEN 'best'
                 ELSE 'inferior'
               END as match_type
        from   alignment a,
               (select query,
                       Max(identity) as top_identity,
                       count(distinct taxon) as num_taxa
                from   alignment
                group  by query) s
        where  a.query = s.query
        group by a.query, a.taxon
        )
group  by taxon,
          match_type; 
" | perl -E '
my %result;
while(<>){
  chomp;
  my ($taxon, $match_type, $c) = split "\t";
  $result{$taxon}{$match_type} = $c}
say join "\t", qw/taxon num_unique num_best num_inferior/;
for my $taxon (sort keys %result){
  say join( "\t", $taxon, (map {$result{$taxon}{$_} // 0} qw/unique best inferior/))
}'

