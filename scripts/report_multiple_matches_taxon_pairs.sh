#!/bin/bash
# Reports multiple matches split by taxon pairs
# Prints pairs of taxa that match in the same query, with better / equal / worse identity
if [ ! "$1" ] ; then
  echo "Usage: $0 database.sql"
  exit 1
fi

sqlite3 -separator $'\t' "$1" <<< "
 select at,
       bt,
       rel,
       Count(*)
from   (select a.taxon at,
               b.taxon bt,
               CASE
                 WHEN a.identity - b.identity between -1e-6 and 1e-6 THEN '=='
                 WHEN a.identity > b.identity THEN '>'
                 ELSE '<'
               END     as rel
        from   alignment a,
               alignment b
        where  a.query = b.query
        and a.taxon != b.taxon
        group by a.query, at, bt, rel
               )
group  by at,
          bt,
          rel
order  by at;
" | perl -E '
my %result;
while(<>){
chomp;
  my ($at, $bt, $rel, $c) = split "\t";
  my $k = "$at\t$bt";
  $result{$k}{$rel} = $c;
}
say join "\t", qw/taxon_1 taxon_2 > == </;
for my $k (sort keys %result){
  say join "\t", $k, map {$result{$k}{$_} // 0} qw/> == </;
}' 
