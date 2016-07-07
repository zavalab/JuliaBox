# Julia script to define some common functions
# Kibaek Kim - ANL/MCS 09252015

# read file and create dictionary
function readDict(f)
	d = readdlm(f);
	return Dict(zip(d[:,1], d[:,2]));
end

# read file and create dictionary in reverse column indices
function readRDict(f)
	d = readdlm(f);
	return Dict(zip(d[:,2], d[:,1]));
end

# create a dictionary for 2-dimensional data
function create2DDict(d)
	dd = Dict{AbstractString, Array{Any,1}}();
	for j in 2:size(d,2)
		setindex!(dd, d[2:25,j], d[1,j]);
	end
	return dd;
end