function tree = j_maketree(I,names)

dist_matrix= 1 ./ 2.^I

tree= seqlinkage(dist_matrix,'UPGMA',names)

