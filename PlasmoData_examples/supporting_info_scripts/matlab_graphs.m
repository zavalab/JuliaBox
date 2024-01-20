rng(10)
dim3 = 1
rand_array = rand(100, 100, dim3)


G = graph()
fmttr = "(%i, %i)"
for i=1:100
    for j=1:100
        G = addnode(G, sprintf(fmttr, i, j))
    end
end


for i=1:100
    for j=1:100
        if (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i, j + 1))
        end
        
        if (i ~= 100)
            G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i + 1, j))
        end
        
        if (i ~= 100) && (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i + 1, j + 1))
        end
        
        if (i ~= 1) && (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i - 1, j + 1))
        end
    end
end

fmttr2 = "weight%i"
for i=1:dim3
    mat_slice = rand_array(:, :, i)
    vec_slice = mat_slice(:)
    G.Nodes.(sprintf(fmttr2, i)) = vec_slice
end

size_val1 = GetSize(G)


rng(10)

dim3 = 10
rand_array = rand(100, 100, dim3)


G = graph()
fmttr = "(%i, %i)"
for i=1:100
    for j=1:100
        G = addnode(G, sprintf(fmttr, i, j))
    end
end


for i=1:100
    for j=1:100
        if (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i, j + 1))
        end
        
        if (i ~= 100)
            G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i + 1, j))
        end
        
        if (i ~= 100) && (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i + 1, j + 1))
        end
        
        if (i ~= 1) && (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i - 1, j + 1))
        end
    end
end

fmttr2 = "weight%i"
for i=1:dim3
    mat_slice = rand_array(:, :, i)
    vec_slice = mat_slice(:)
    G.Nodes.(sprintf(fmttr2, i)) = vec_slice
end

size_val10 = GetSize(G)

rng(10)

dim3 = 25
rand_array = rand(100, 100, dim3)


G = graph()
fmttr = "(%i, %i)"
for i=1:100
    for j=1:100
        G = addnode(G, sprintf(fmttr, i, j))
    end
end


for i=1:100
    for j=1:100
        if (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i, j + 1))
        end
        
        if (i ~= 100)
            G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i + 1, j))
        end
        
        if (i ~= 100) && (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i + 1, j + 1))
        end
        
        if (i ~= 1) && (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i - 1, j + 1))
        end
    end
end

fmttr2 = "weight%i"
for i=1:dim3
    mat_slice = rand_array(:, :, i)
    vec_slice = mat_slice(:)
    G.Nodes.(sprintf(fmttr2, i)) = vec_slice
end

size_val25 = GetSize(G)

rng(10)

dim3 = 100
rand_array = rand(100, 100, dim3)


G = graph()
fmttr = "(%i, %i)"
for i=1:100
    for j=1:100
        G = addnode(G, sprintf(fmttr, i, j))
    end
end


for i=1:100
    for j=1:100
        if (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i, j + 1))
        end
        
        if (i ~= 100)
            G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i + 1, j))
        end
        
        if (i ~= 100) && (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i + 1, j + 1))
        end
        
        if (i ~= 1) && (j ~= 100)
           G = addedge(G, sprintf(fmttr, i, j), sprintf(fmttr, i - 1, j + 1))
        end
    end
end

fmttr2 = "weight%i"
for i=1:dim3
    mat_slice = rand_array(:, :, i)
    vec_slice = mat_slice(:)
    G.Nodes.(sprintf(fmttr2, i)) = vec_slice
end

size_val100 = GetSize(G)

function size = GetSize(this) 
   props = properties(this); 
   totSize = 0; 
   
   for ii=1:length(props) 
      currentProperty = getfield(this, char(props(ii))); 
      s = whos('currentProperty'); 
      totSize = totSize + s.bytes; 
   end
  
   %fprintf(1, '%d bytes\n', totSize); 
   size = totSize
end

