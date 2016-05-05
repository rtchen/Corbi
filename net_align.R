net_align <- function(netfileA, netfileB,simfile, alpha=1,beta=2,delta.d=1e-10,output="result.txt")
{
  combined_net <- read_net(netfileA,netfileB,simfile)
  crf <-build_model(combined_net,node_sim,alpha,beta,delta.d)
  result <- decode.lbp(crf)
  result <- crf$state.map[cbind(1:crf$n.nodes, result)]
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
  net.node <- c(net.nodeA,net.nodeB)
  net.size = length(net.node)
  net.edgeA <- cbind(as.character(net.textA[,1]), as.character(net.textA[,2]))
  net.edgeB <- cbind(as.character(net.textB[,1]), as.character(net.textB[,2]))
  net.edgeA <- net.edgeA[net.edgeA[,2] != "", ]
  net.edgeB <- net.edgeB[net.edgeB[,2] != "", ]
  sim.edge <- cbind(as.character(sim.text[,1]), as.character(sim.text[,2]))
  net.edge <- rbind(net.edgeA,sim.edge,net.edgeB)
  node.id <- seq_along(net.node)
  names(node.id) <- net.node
  net.matrix <- sparseMatrix(node.id[net.edge[,1]], node.id[net.edge[,2]], x=T, dims=c(net.size, net.size))
  net.matrixA <- sparseMatrix(node.id[net.edgeA[,1]], node.id[net.edgeA[,2]], x=T, dims=c(net.size+1, net.size+1))
  net.matrixB <- sparseMatrix(node.id[net.edgeB[,1]], node.id[net.edgeB[,2]], x=T, dims=c(net.size+1, net.size+1))
  sim.matrix <- sparseMatrix(node.id[sim.edge[,1]], node.id[sim.edge[,2]], x=as.numeric(sim.text[,3]), dims=c(net.size, net.size))
  list(node_sim = sim.matrix,node.id = node.id,sizeA=net.sizeA,node=net.node, matrix=net.matrix,
       matrixA = net.matrixA,matrixB = net.matrixB)
  
}

build_model <- function(combined_net,node_sim,alpha,beta,delta.d)
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
    W <- as.matrix(W2[crf$state.map[n1, m1], crf$state.map[n2, m2]])
    crf$edge.pot[[e]] <- exp(W*beta/2)
    }
    else if(n1 > A.size && n2 > A.size){
    W <- as.matrix(W1[crf$state.map[n1, m1], crf$state.map[n2, m2]])
    crf$edge.pot[[e]] <- exp(W*beta/2)
    }
    else{
      G <- matrix(0,nrow=length(m1),ncol = length(m2))
      for(m11 in m1){
        for(m22 in m2){
          if(crf$state.map[n1,m11]==n2 && crf$state.map[n2,m22]==n1){
            G[m11,m22]<-1
          }
          else if(crf$state.map[n1,m11]!=n2 && crf$state.map[n2,m22]!=n1){
            G[m11,m22]<-1
          }
        }
      }
      crf$edge.pot[[e]]<-G
    }
  }
  crf
}



write_result <- function(combined_net, result, filename="result.txt")
{   
  con <- file(as.character(filename), "w")
  label.name <- c(combined_net$node, "gap")
  for (i in 1:combined_net$sizeA)
  { 
    writeLines(paste(combined_net$node[i], "  ", label.name[result[i]]), con, sep="\n")
  }
  writeLines("", con, sep="\n")
  close(con)
}


