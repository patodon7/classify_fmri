
##############################################################
######Functional Network Creation
##############################################################
rm(list = ls())

##################Change these lines##########################
setwd("/path/to/my/directory")
####Read in list of filenames- run ICA/make distance matrices
filenames<-read.table("myfilenames.txt")
##############################################################
num_subjects<-dim(filenames)[1]





####Loop through to run ICA on each scan and create a distance matrix calculating the temporal correlation (or funtional connectivity) between pairs of components.  
####Call the library used for analysis.  This must be preinstalled.
library(AnalyzeFMRI)
###The number of components is set as default to 20.
num_components<-20
###Now create the array that will hold the "distance" matrix for each subject.    
data_array_distance<-array(NA, c(num_subjects, num_components, num_components))
for (i in 1:num_subjects)
{	#Perform ICA on each file
	temp<-f.ica.fmri(as.character(filenames[i,1]), n.comp = num_components)
	for(j in 1:num_components)
		{	for(k in j:num_components)
			{	#Calculate the "distance" between two components as the maximal cross-correlation taken over a timelag.
				data_array_distance[i,k,j]<-max(1-abs(max(ccf(temp$A[,j], temp$A[,k], plot = FALSE, lag.max = 10)$acf,0)))
				data_array_distance[i,j,k]<-data_array_distance[i,k,j]	}	}	
	diag(data_array_distance[i,,])<-0}
s<-rep(num_components, num_subjects)


###If using melodic_mix files instead:

##################Change these lines##########################
setwd("/path/to/my/directory")
filenames<-read.table("files_melodic_mix.txt")
##############################################################

###Find max numbers of components
s<-c(rep(0,num_subjects))
for (i in 1:num_subjects)
{	s[i]<-dim((read.table(as.character(filenames[i,1]))))[2]	}
###Now read in the data and create the distance matrix.
data_array_distance<-array(NA, c(num_subjects, max(s), max(s)))
for (i in 1:num_subjects)
{	##Read in ICA results
	temp<-as.matrix(read.table(as.character(filenames[i,1])))
	for(j in 1:s[i])
	{	for(k in j:s[i])
		{	#Calculate the "distance" between two components as the maximal cross-correlation taken over a timelag.
			data_array_distance[i,k,j]<-1-abs(max(ccf(temp[,j], temp[,k], plot = FALSE, lag.max = 10)$acf))
			data_array_distance[i,j,k]<-data_array_distance[i,k,j]	}	}	
	diag(data_array_distance[i,,])<-0		}


##############################################################
######Feature Extraction
##############################################################

###Turn each distance matrix into a graph, and extract summary properties of the graph
library(igraph)
library(vegan)

####Function to create "graph" object.
makegraph<-function(my_iso)
{	##dim is dimension of matrix
	my_dist<-as.matrix(dist(my_iso$points[]))
	k<-dim(my_dist)[1]
	my_net<-matrix(0, nrow = k, ncol = k)
	which.rows<-my_iso$net[,1]
	which.cols<-my_iso$net[,2]
	for(j in 1:length(which.rows))
	{	my_net[which.rows[j], which.cols[j]] <-my_dist[which.rows[j], which.cols[j]]
	my_net[which.cols[j], which.rows[j]] <-my_dist[which.cols[j], which.rows[j]]		}	
	my_net
}



my_feature_matrix<-matrix(NA, nrow = num_subjects, ncol = 12)
for(i in 1:num_subjects)
{	d<-matrix(data_array_distance[i,1:s[i],1:s[i]], nrow = s[i])	
	my_iso<-isomap(d[1:s[i],1:s[i]],axes=1, epsilon=median(d), ndim = 10)
	my_net<-makegraph(my_iso)
	d2 <- graph.adjacency(my_net, weighted = TRUE )	
	my_feature_matrix[i,]<-c(average.path.length(d2),clique.number(d2),graph.density(d2),edge.connectivity(d2),median(closeness(d2)) ,median(graph.coreness(d2)),max(degree(d2)),median(degree(d2)),min(degree(d2)),vcount(d2),ecount(d2),transitivity(d2))
	}

colnames(my_feature_matrix)<-c("Average Path Length", "Clique Number", "Graphy Density", "Edge Connectivity", "Median Closeness", "Graph Coreness", "Max Degree", "Median Degree", "Min Degree", "Vertex Count", "Edge Count", "Transitivity")


##############################################################
######SVM Classification and Cross-Validation
##############################################################
library(e1071)
my_cat<-c(rep("Normal",17), rep("Schizophrenic",10))
my_cat<-as.factor(my_cat)
my_svm<-svm(my_cat~., data = my_feature_matrix, cross=10)
my_svm$accuracies
