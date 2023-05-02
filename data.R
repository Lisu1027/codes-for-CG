## Function for constructing diagonal precision matrix
precision_matrix = function(c1=0.9,c2=0.6,c3=0,p1,p2,p3){
  p = sum(p1)+sum(p2)+sum(p3)
  matrix = matrix(0,p,p)
  for(i in (1:length(p1))){
    if(i == 1){
      n1 = 0
      n2 = p1[i]
      matrix[1:n2,1:n2] = matrix(c1,n2,n2)
    }
    if(i != 1){
      n1 = sum(p1[1:(i-1)])
      n2 = p1[i]
      matrix[(n1+1):(n1+n2),(n1+1):(n1+n2)] = matrix(c1,n2,n2)
    }
  }
  for(i in (1:length(p2))){
    if(i == 1){
      n1 = sum(p1)
      n2 = p2[i]
      matrix[(n1+1):(n1+n2),(n1+1):(n1+n2)] = matrix(c2,n2,n2)
    }
    if(i != 1){
      n1 = sum(p1)+sum(p2[1:(i-1)])
      n2 = p2[i]
      matrix[(n1+1):(n1+n2),(n1+1):(n1+n2)] = matrix(c2,n2,n2)
    }
  }
  for(i in (1:length(p3))){
    if(i == 1){
      n1 = sum(p1)+sum(p2)
      n2 = p3[i]
      matrix[(n1+1):(n1+n2),(n1+1):(n1+n2)] = matrix(c3,n2,n2)
    }
    if(i != 1){
      n1 = sum(p1)++sum(p2)+sum(p3[1:(i-1)])
      n2 = p3[i]
      matrix[(n1+1):(n1+n2),(n1+1):(n1+n2)] = matrix(c3,n2,n2)
    }
  }
  diag(matrix) = 1
  return(matrix)
}

## Function for constructing a simple off-diagonal precision matrix
precision_matrix_dno = function(precision_matrix_d,c,p1){
  for(j in (2:(length(p1)-1))){
    if(j == 2){a11=1}
    if(j != 2){a11 = sum(p1[1:(j-2)])+1}
    a12 = sum(p1[1:(j-1)])
    b11 = sum(p1[1:j])+1
    b12 = sum(p1[1:(j+1)])
    precision_matrix_d[a11:a12,b11:b12] = c
    precision_matrix_d[b11:b12,a11:a12] = c
  }
  return(precision_matrix_d)
}
