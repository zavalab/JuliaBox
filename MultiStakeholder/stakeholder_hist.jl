# Plot histogram of stakeholder dissatisfactions
# Use with three_objective2.jl
# Created by Alex Dowling (adowling@wisc.edu)
# at the University of Wisconsin-Madison
# Last Modified: Feb. 17th, 2016

figure()

# Font size for annotation
fsa = 18

# Font size for alpha
fsal = 18

# Font size for axis labels
fsax = 18

# Max size for y
yMax = 120

# Height for alpha
alphaH = 85

# Height for median label
medianH = 65

subplot(311)
dis = sv1[1].z - zp1
plt[:hist](dis,35,alpha=0.75,color="red")
xlim(0,1)
ylim(0,yMax)

annotate("Largest", fontsize=fsa, xy=[maximum(dis),0], xycoords="data",
            xytext=(maximum(dis)-0.1, 60), textcoords="data",
            arrowprops=Dict("facecolor"=>"black", "shrink"=>0.05),
            horizontalalignment="right", verticalalignment="top",
            )

md = median(dis)
plot([md,md],[0,150],color="black",linewidth=3,linestyle="--")
text(md+0.02,medianH,"Median",fontsize=fsa)
            
text(0.9, alphaH,L"\alpha = 0",va="center",ha="center",fontsize=fsal,color="red")

subplot(312)
dis = sv1[51].z - zp1
plt[:hist](dis,25,alpha=0.75,color="green")
xlim(0,1)
ylim(0,yMax)
ylabel("Number of Stakeholders",fontsize=fsax)

annotate("Largest", fontsize=fsa, xy=[maximum(dis),0], xycoords="data",
            xytext=(maximum(dis)+0.1, 60), textcoords="data",
            arrowprops=Dict("facecolor"=>"black", "shrink"=>0.05),
            horizontalalignment="left", verticalalignment="top",
            )
text(0.9, alphaH,L"\alpha = \frac{1}{2}",va="center",ha="center",fontsize=fsal,color="green")

md = median(dis)
plot([md,md],[0,150],color="black",linewidth=3,linestyle="--")
text(md+0.02,medianH,"Median",fontsize=fsa)

subplot(313)
dis = sv1[101].z - zp1
plt[:hist](dis,25,alpha=0.75,color="blue")
xlim(0,1)
ylim(0,yMax)

annotate("Largest", fontsize=fsa, xy=[maximum(dis),0], xycoords="data",
            xytext=(maximum(dis)+0.1, 60), textcoords="data",
            arrowprops=Dict("facecolor"=>"black", "shrink"=>0.05),
            horizontalalignment="left", verticalalignment="top",
            )
            
text(0.9, alphaH,L"\alpha = 1",va="center",ha="center",fontsize=fsal,color="blue")

md = median(dis)
plot([md,md],[0,150],color="black",linewidth=3,linestyle="--")
text(md+0.02,medianH,"Median",fontsize=fsa)

xlabel("Stakeholder Dissatisfactions",fontsize=fsax)

# Why is the largest dissatisfaction for alpha = 0.5 smaller than the largest dissatisfaction for alpha = 1.0?