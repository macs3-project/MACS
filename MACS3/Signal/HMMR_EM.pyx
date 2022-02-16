# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-02-16 13:58:11 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
import numpy as np
#from scipy.stats import norm    # there is another implemented function MACS3.Signal.Prob.pnorm
from MACS3.Signal.Prob import pnorm

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
cdef inline float get_weighted_density( x, mean, stddev, weights ):
    """Description:
    
    parameters:
    return value:
    """
    return weight * pnorm( x, mean, stddev )
    # return np.multiply(weight, normal_dist_density) # make sure both are np.array types

cdef return_greater( self, list data ):
    """
    Return the index of the largest value in an array of doubles
    @param data: an Array of doubles 
    @return an integer representing thre index of the largest value in the inputted array
    """    
    # TL: there must be a easy function in numpy/scipy or even Python itself for this
    # also, this function doesn't use any 'self', so could be a general function outside of this class
    largest_index = -1
    greatest_value = -1.0
    for i in range(0, len(data)):
        if data[i] > greatest_value:
            greatest_value = data[i]
            largest_index = i
    for i in range(0, len(data)): # if there are other mostly likely category, ignore this data point
        if i != largest_index:
            if data[i] == greatest_value:
                largest_index = -1
    return largest_index
    
cdef class HMMR_EM:
    """ Main HMMR EM class.
    
    """
    cdef:
        # define the values that should be returned
        list fragMeans
        list fragStddevs
        # private variables
        float episilon = 0.0005
        int maxIter = 20
        float jump = 1.5
        object bamfile # a MACS3.IO.BAM.BAMAccessor object

        bool converged = False # when True, terminate the iterations
        np.ndarray[np.int32_t, ndim=1] data # data for fragment lengths
        list mu
        list lamba
        list weight
        
    def __init__ ( self, object petrack, list em_means, list em_stddevs , dict genome, float sample_percentage ):
        """Initialize HMMR_EM object.

        parameters:
            1. petrack: a MACS3.Signal.PairedEndTrack.PETrackI object
            1. em_means: list of initial means of fragments, for short, mono, di, and tri signals
            2. em_stddevs: list of initial stddevs of fragments, for short, mono, di, and tri signals
            3. genome: a dictionary containing the lengths (as value) of each chromosomes (as key)
        """
    
        # step1: initialize the values of fragMeans and fragStddevs
        # * Set fragment length distribution parameters. 
        #  * Use inputed values to set initial values, if provided. 
        #  * Else use defaults
        #  */
        # initial values
        self.petrack = petrack # we may need to use a deepcopy
        self.epsilon = 0.0005
        self.maxIter = 20
        self.jump = 1.5
        self.fragMeans = np.array(em_means)
        self.fragStddevs = np.array(em_stddev)
        self.sample_percentage = sample_percentage

        # first, let's prepare the lengths data
        # sample down
        self.petrack.sample_percent( self.sample_percentage, seed=12345 ) # may need to provide seed option for init function
        self.data = self.petrack.fraglengths()
        # then we only keep those with fragment lengths > 100 and < 1000
        self.data = self.data[ self.data > 100 ]
        self.data = self.data[ self.data < 1000 ]
        
        # next, we will calculate the weights -- ie the proportion of fragments in each length category
        cutoff1 = em_means[ 1 ] - em_means[ 0 ]/2 + em_means[ 0 ]
        cutoff2 = em_means[ 2 ] - em_means[ 1 ]/2 + em_means[ 1 ]
        cutoff3 = em_means[ 3 ] - em_means[ 2 ]/2 + em_means[ 2 ]
        
        sum4 = len( self.data )
        sum3 = sumn( self.data < cutoff3 )
        sum2 = sumn( self.data < cutoff2 )
        sum1 = sumn( self.data < cutoff1 )

        counter4 = sum4 - sum3
        counter3 = sum3 - sum2
        counter2 = sum2 - sum1
        counter1 = sum1

        self.weight = [ counter1/sum4, counter2/sum4, counter3/sum4, counter4/sum4]

        self.__learn()
        #self.__getMeans()
        #self.__getLambda()
	# double[] tempMeans = em.getMeans();
	# double[] tempLam = em.getLamda();
	# em = null;
	# for (int i = 0;i < tempMeans.length;i++){
	# 	//This will update the parameters IFF they were updated. If they become NaN, leave as default
	# 	if(!Double.isNaN(tempMeans[i]) && !Double.isNaN(tempLam[i])){
	# 		fragMeans[i+1] = tempMeans[i];
	# 		fragStddevs[i+1] = tempLam[i];
	# 	}
	# }
	# em = null;lengths = null;tempMeans=null;tempLam=null;

        # in the end, we return. And users should be able to retrieve
        # fragMeans and fragStddevs.
        return


    #cdef void __pullLargeLengths(self, object bamfile, options.min_map_quality, GENOME, options.em_means):
        # TL: I will try to implement this
# public class pullLargeLengths {
	
# 	private File bam;
# 	private File index;
# 	private int minQ;
# 	private double[] lengths;
# 	private ArrayList<TagNode> genome;
# 	private double[] weights = new double[3];
# 	private double[] means;
# 	/**
# 	 * Constructor for creating pullLargeLengths object, reading the data and setting weights
# 	 * @param b a File representing the BAM file of reads
# 	 * @param i a File representing the BAM index file of reads
# 	 * @param m an integer representing the minimum mapping quality of reads to pull
# 	 * @param g an ArrayList of TagNode representing the genome
# 	 * @param mu an Array of doubles representing the means of the fragment length distributions
# 	 */
# 	public pullLargeLengths(File b, File i,int m,ArrayList<TagNode> g,double[] mu){
# 		bam = b;
# 		index = i;
# 		minQ = m;
# 		genome = g;
# 		means = mu;
# 		read();
# 		setWeights();
# 	}
# 	/**
# 	 * Access the lengths
# 	 * @return an Array of doubles representing the lengths of the fragments read
# 	 */
# 	public double[] getLengths(){return lengths;}
# 	/**
# 	 * Access the weights
# 	 * @return an Array of doubles representing the proportion of fragments in each length category
# 	 */
# 	public double[] getWeights(){return weights;}
# 	/**
# 	 * Access the down sampled weights
# 	 * @param div an integer representing the factor to divide the number of lengths by to downsample
# 	 * @return an Array of doubles representing the downsampled lengths
# 	 */
# 	public double[] getSampledLengths(int div){
# 		int len = lengths.length/div;
# 		double[] newLengths = new double[len];
# 		shuffle(lengths);
# 		for(int i = 0; i < len;i++){
# 			newLengths[i] = lengths[i];
# 		}
# 		return newLengths;
# 	}
# 	/**
# 	 * Shuffle the length data
# 	 * @param list an Array of doubles representing the pulled fragment lengths
# 	 */
# 	private void shuffle(double[] list){
# 		Random rnd = new Random();
# 		for(int i = list.length -1;i > 0;i--){
# 			int index = rnd.nextInt(i+1);
# 			double a = list[index];
# 			list[index] = list[i];
# 			list[i] = a;
# 		}
# 	}
# 	/**
# 	 * Read the data and create a list of lengths
# 	 */
# 	private void read(){
# 		int counter = 0;
# 		SAMFileReader reader = new SAMFileReader(bam,index);
# 		ArrayList<Double> temp = new ArrayList<Double>();
# 		for (int i = 0; i < genome.size();i++){
# 			String chr = genome.get(i).getChrom();
# 			int start = genome.get(i).getStart();
# 			int stop = genome.get(i).getStop();
# 			CloseableIterator<SAMRecord> iter = reader.query(chr,start,stop,false);
# 			while (iter.hasNext()){
# 				SAMRecord record = null;
# 				try{
# 					record = iter.next();
# 				}
# 				catch(SAMFormatException ex){
# 					System.out.println("SAM Record is problematic. Has mapQ != 0 for unmapped read. Will continue anyway");
# 				}
# 				if(record != null){
# 				if(!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && record.getFirstOfPairFlag()) {
# 					if (record.getMappingQuality() >= minQ){
						
# 						if (Math.abs(record.getInferredInsertSize()) > 100 && Math.abs(record.getInferredInsertSize())
# 								< 1000){
# 							counter+=1;
# 							temp.add((double)Math.abs(record.getInferredInsertSize()));
# 						}
# 					}
# 				}
# 				}
# 			}
# 			iter.close();
# 		}
# 		reader.close();
# 		lengths = new double[counter];
# 		for (int i = 0;i < temp.size();i++){
# 			if (temp.get(i) > 100){
# 				lengths[i] = temp.get(i);
# 			}
# 		}
		
# 	}
# 	/**
# 	 * Set the weights ie the proportion of fragments in each length category
# 	 */
# 	private void setWeights(){
# 		double cutone = (means[1] - means[0])/2 + means[0];
# 		double cuttwo = (means[2] - means[1])/2 + means[1];
# 		int counter1=0;
# 		int counter2=0;
# 		int counter3=0;
# 		for (int i = 0;i < lengths.length;i++){
# 			if (lengths[i] < cutone)
# 				counter1++;

    cdef __iterate(self):
        """Description: This is a private function only used by HMMR_EM class

        parameters:
        return value:
        """
        cdef:
            list temp, counter, means, stds
            float total
            int i, j
            
        temp = []               # for each category, the likelihood
        counter = []            # for each category, the number of data points/fragment
        total = 0.0             # total number of data points/fragments assigned to four categories
        means = []              # for each category, the new mean
        stds = []               # for each category, the new stddev
        for i in range( 0, len( self.data ) ):
            for j in range( 0, len( self.mu ) ):
                # for each category: short, mono, di, tri- (4 in total)
                temp[j] = get_weighted_density( self.data[i], self.mu[j], self.lam[j], self.weights[j] )
            # now look for the most likely category, as `index`
            index = return_greater( temp )
            
            # is this too large of a file to save means,stds for each iteration?
            if index != -1: # If we can find a mostly likely category
                ##----
                means[index] = np.mean( self.data[0:i] ) #check - do we want means[index] or means[i] or means.append()
                stds[index] = np.std(self.data[0:i])
                ##---- This block is wrong. We need to implement incremental way to calculate mean and stddev
                ## https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
                ## Check this
                counter[index] += 1.0 # check on this - do we want it to be a list with increasing values?
                total += 1.0
        for i in range( 0, len( self.mu ) ): # what are the lengths of means, stds, data, mu?
            # we will emplify the difference between new and old means/stds, and update weights.
            self.mu[i] = self.mu[i] + (jump * (means[i] - self.mu[i]))
            self.lam[i] = self.lam[i] + (jump * (stds[i] - self.lam[i]))
            self.weights[i] = counter[i] / total

    cdef bool learn(self):
        """Description:

        parameters:
        return value:
        """
        cdef:
            int itr = 0         # number of iterations
            int i
            int counter         # number of category that has been converged
        self.converged = False
        while self.converged == False:
            for i in range( 0, len( self.mu ) ):
                old_mu[a] = self.mu[a]
                old_lam[a] = self.lam[a]
                old_weights[a] = self.weights[a]

            self.__iterate()

            counter = 0
            for i in range(0, len(self.mu)):
                if abs(old_mu[a] - self.mu[a]) < self.epsilon) and abs(old_weights[a] - self.weights[a]) < self.epsilon) and abs(old_lam[a] - self.lam[a]) < self.epsilon)
                    counter += 1
            if counter == len( self.mu ):
                self.converged = True
            itr += 1
            if itr >= self.maxIter:
                break
        return self.converged
            
