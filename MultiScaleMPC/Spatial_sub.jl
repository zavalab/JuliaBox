function SetIndexSet(gridSize)
  totalIndex=Set{Tuple{Int64,Int64}}()
  neighborIndex=Set{Tuple{Int64,Int64}}()
  innerIndex=Set{Tuple{Int64,Int64}}()

  for i=1:gridSize
    push!(neighborIndex,(0,i))
    push!(neighborIndex,(gridSize+1,i))
    push!(neighborIndex,(i,0))
    push!(neighborIndex,(i,gridSize+1))
  end
  for i=1:gridSize
    for j=1:gridSize
      push!(innerIndex,(i,j))
    end
  end
  totalIndex=union(innerIndex,neighborIndex)
  neighbor=Dict(ind=>Set([(ind[1]-1,ind[2]),(ind[1]+1,ind[2]),(ind[1],ind[2]-1),(ind[1],ind[2]+1)]) for ind in innerIndex)

  return totalIndex,neighborIndex,innerIndex,neighbor
end

function SetPartitionIndexSet(partitionPerSide,partitionSize)

  totalPartitionIndex=Set{Tuple{Int64,Int64}}()
  for p=1:partitionPerSide
    for q=1:partitionPerSide
      push!(totalPartitionIndex,(p,q))
    end
  end
  totalIndexPart=Dict(part=>Set{Tuple{Int64,Int64}}() for part in totalPartitionIndex)
  neighborIndexPart=Dict(part=>Set{Tuple{Int64,Int64}}() for part in totalPartitionIndex)
  boundaryIndexPart=Dict(part=>Set{Tuple{Int64,Int64}}() for part in totalPartitionIndex)
  innerIndexPart=Dict(part=>Set{Tuple{Int64,Int64}}() for part in totalPartitionIndex)
  partitionTransfer=Dict(ind=>(0,0) for ind in totalIndex)

  for part in totalPartitionIndex
    p=part[1]
    q=part[2]
    X=partitionSize*(p-1)
    Y=partitionSize*(q-1)

    for i=1:partitionSize
      push!(neighborIndexPart[part],(X,Y+i))
      push!(neighborIndexPart[part],(X+partitionSize+1,Y+i))
      push!(neighborIndexPart[part],(X+i,Y))
      push!(neighborIndexPart[part],(X+i,Y+partitionSize+1))

      push!(boundaryIndexPart[part],(X+1,Y+i))
      push!(boundaryIndexPart[part],(X+partitionSize,Y+i))
      push!(boundaryIndexPart[part],(X+i,Y+1))
      push!(boundaryIndexPart[part],(X+i,Y+partitionSize))
    end
    for i=1:partitionSize
      for j=1:partitionSize
        push!(innerIndexPart[part],(X+i,Y+j))
        partitionTransfer[(X+i,Y+j)]=part
      end
    end
    totalIndexPart[part]=union(innerIndexPart[part],neighborIndexPart[part])
  end



  return totalPartitionIndex,totalIndexPart,neighborIndexPart,boundaryIndexPart,innerIndexPart,partitionTransfer
end

function SetCoarseIndexSet(coarseSize,partitionSize,totalPartitionIndex,totalIndexPart)
  coarsePerSide=partitionSize/coarseSize
  coarseIndexPart=Dict{Tuple{Int64,Int64},Dict{Tuple{Int64,Int64},Set{Tuple{Int64,Int64}}}}()
  coarseInnerIndex=Set{Tuple{Int64,Int64}}()
  coarseNeighborIndex=Set{Tuple{Int64,Int64}}()
  coarseTransfer=Dict(part=>Dict{Tuple{Int64,Int64},Tuple{Int64,Int64}}() for part in totalPartitionIndex)
  for r=1:coarsePerSide
    push!(coarseNeighborIndex,(0,r))
    push!(coarseNeighborIndex,(coarsePerSide+1,r))
    push!(coarseNeighborIndex,(r,0))
    push!(coarseNeighborIndex,(r,coarsePerSide+1))
  end
  for r=1:coarsePerSide
    for s=1:coarsePerSide
      push!(coarseInnerIndex,(r,s))
    end
  end
  coarseTotalIndex=union(coarseNeighborIndex,coarseInnerIndex)

  for part in totalPartitionIndex
    coarseIndexPart[part]=Dict(coarse=>Set{Tuple{Int64,Int64}}() for coarse in coarseTotalIndex)
    coarseTransfer[part]=Dict(ind=>(0,0) for ind in totalIndexPart[part])
    for ind in totalIndexPart[part]
      r=div(ind[1]-(part[1]-1)*partitionSize-1+coarseSize,coarseSize)
      s=div(ind[2]-(part[2]-1)*partitionSize-1+coarseSize,coarseSize)
      push!(coarseIndexPart[part][(r,s)],ind)
      coarseTransfer[part][ind]=(r,s)
    end
  end
  return coarseIndexPart,coarseTotalIndex,coarseInnerIndex,coarseNeighborIndex,coarseTransfer
end

function SetOrder1()
  order=[]
  for p=1:partitionPerSide
    if mod(p,2)==1
      for q=1:partitionPerSide
        push!(order,(p,q))
      end
    end
    if mod(p,2)==0
      for q=partitionPerSide:-1:1
        push!(order,(p,q))
      end
    end
  end
  return order
end

function SetOrder4()
  dPartition=Dict(part=>0.0 for part in totalPartitionIndex)
  for part in totalPartitionIndex
    dSum=0
    for ind in innerIndexPart[part]
      dSum=dSum+d[ind]^2
    end
    dPartition[part]=dSum
  end
  preOrder=sort(collect(dPartition), by=x->x[2],rev=true)
  order=Array{Tuple{Int64,Int64}}(totNumPartitions)
  for i=1:totNumPartitions
    order[i]=preOrder[i][1]
  end
  return order
end



function SetOrder3()
  order=[]
  for p=1:partitionPerSide
    for q=1:partitionPerSide
      if mod(p+q,2)==0
        push!(order,(p,q))
      end
    end
  end
  for p=1:partitionPerSide
    for q=1:partitionPerSide
      if mod(p+q,2)==1
        push!(order,(p,q))
      end
    end
  end
  return order
end

function SetOrder2()
  p=1
  q=1
  order=[]
  remains=copy(totalPartitionIndex)
  push!(order,(1,1))
  pop!(remains,(1,1))
  while length(remains)!=0
    while in((p+1,q),remains)
      p+=1
      push!(order,(p,q))
      pop!(remains,(p,q))
    end
    while in((p,q+1),remains)
      q+=1
      push!(order,(p,q))
      pop!(remains,(p,q))
    end
    while in((p-1,q),remains)
      p-=1
      push!(order,(p,q))
      pop!(remains,(p,q))
    end
    while in((p,q-1),remains)
      q-=1
      push!(order,(p,q))
      pop!(remains,(p,q))
    end
  end
  return order
end
