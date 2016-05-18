net_align <- function(netfileA, netfileB,simfile, alpha=1,beta=2,delta.d=1e-10,output="result.txt")
{ a <- Sys.time()
  combined_net <- read_net(netfileA,netfileB,simfile)
  b<- Sys.time()
  crf <-build_model(combined_net,alpha,beta,delta.d)
  c <- Sys.time()
  result <- decode.lbp(crf)
  d <- Sys.time()
  result <- crf$state.map[cbind(1:crf$n.nodes, result)]
  print(b-a)
  print(c-b)
  print(d-c)
  write_result(combined_net, result, output)
}


read_net <- function(netfileA,netfileB,simfile)
{ 
  net.textA <- as.matrix(read.table(netfileA, fill=T, as.is=T))
  net.textB <- as.matrix(read.table(netfileB, fill=T, as.is=T))
  sim.text <- as.matrix(read.table(simfile, fill=T, as.is=T))
  net.nodeA <- unique(as.character(net.textA))
  net.nodeA <- net.nodeA[net.nodeA != ""]
  net.sizeA <- length(net.nodeA)
  net.nodeB <- unique(as.character(net.textB))
  net.sizeB <- length(net.nodeB)
  net.nodeB <- net.nodeB[net.nodeB != ""]
  sim.node <- unique(as.character(sim.text))
  net.node <- unique(c(net.nodeA,net.nodeB,sim.text))
  net.size = length(net.node)
  net.edgeA <- cbind(as.character(net.textA[,1]), as.character(net.textA[,2]))
  net.edgeB <- cbind(as.character(net.textB[,1]), as.character(net.textB[,2]))
  net.edgeA <- net.edgeA[net.edgeA[,2] != "", ]
  net.edgeB <- net.edgeB[net.edgeB[,2] != "", ]
  sim.edge <- cbind(as.character(sim.text[,1]), as.character(sim.text[,2]))
  net.edge <- rbind(net.edgeA,sim.edge,net.edgeB)
  node.id <- seq_along(net.node)
  names(node.id) <- net.node
  net.matrixA <- sparseMatrix(node.id[net.edgeA[,1]], node.id[net.edgeA[,2]], x=T, dims=c(net.size+1, net.size+1))
  net.matrixB <- sparseMatrix(node.id[net.edgeB[,1]], node.id[net.edgeB[,2]], x=T, dims=c(net.size+1, net.size+1))
  sim.matrix <- sparseMatrix(node.id[sim.edge[,1]], node.id[sim.edge[,2]], x=as.numeric(sim.text[,3]), dims=c(net.size, net.size))
  net.matrix <- sparseMatrix(node.id[net.edge[,1]], node.id[net.edge[,2]], x=T, dims=c(net.size, net.size))
  list(node_sim = sim.matrix,sizeA=net.sizeA,node=net.node, matrix=net.matrix,size = net.size,
       matrixA = net.matrixA,matrixB = net.matrixB)
  
}

build_model <- function(combined_net,alpha,beta,delta.d)
{ 
  S <- cbind(combined_net$node_sim,delta.d)
  crf <- make.crf(combined_net$matrix, rowSums(S>0)) 
  crf$state.map <- matrix(0, nrow=crf$n.nodes, ncol=crf$max.state)
  for (i in 1:crf$n.nodes)
  { 
    crf$state.map[i, 1:crf$n.states[i]] <- which(S[i,] > 0)
    crf$node.pot[i, 1:crf$n.states[i]] <- exp(S[i, crf$state.map[i,]]*alpha/2)
  }
  A.size = combined_net$sizeA
  W1 <-  combined_net$matrixA
  W2 <- combined_net$matrixB

  
  for (e in 1:crf$n.edges)
  {
    n1 <- crf$edges[e, 1]
    n2 <- crf$edges[e, 2]
    m1 <- 1:crf$n.states[n1]
    m2 <- 1:crf$n.states[n2]
    if(n1<=A.size && n2 <=A.size){
    W <- submatrix(W2,crf$state.map[n1, m1],crf$state.map[n2, m2])
    crf$edge.pot[[e]] <- exp(W*beta/2)
    }
    else if(n1 > A.size && n2 > A.size){
    W <- submatrix(W1,crf$state.map[n1, m1],crf$state.map[n2, m2])
    crf$edge.pot[[e]] <- exp(W*beta/2)
    }
    else{
      G <- matrix(1,nrow=length(m1),ncol = length(m2))
      m11 <- which(crf$state.map[n1,]==n2)
      G[m11,m2]<-0
      m22 <- which(crf$state.map[n2,]==n1)
      G[m1,m22]<-0
      G[m11,m22] <-1
      crf$edge.pot[[e]]<-G
    }
  }
  crf
}



write_result <- function(combined_net, result, filename="result.txt")
{   
  con <- file(as.character(filename), "w")
  label.name <- c(combined_net$node, "gap")
  cnt <- 0
  for (i in 1:combined_net$sizeA)
  { if (result[i] <= combined_net$size){
    cnt <- cnt+1
     }
    writeLines(paste(combined_net$node[i], "  ", label.name[result[i]]), con, sep="\n")
  }
  writeLines("", con, sep="\n")
  print(cnt)
  close(con)
}

column <- function(m, i)
{
   if (inherits(m, "CsparseMatrix")) {
    v <- vector(typeof(m@x), m@Dim[1])
    p <- (m@p[i]+1):m@p[i+1]
    if (p[1] <= p[length(p)])
      v[m@i[p]+1] <- m@x[p]
  }
  else
    v <- m[,i]
  v
}


submatrix <- function(m, rows, cols)
{
  sapply(cols, function(i) column(m, i)[rows])
}
