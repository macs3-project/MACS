# cython: language_level=3
# cython: profile=True
# Time-stamp: <2022-02-02 14:28:52 Tao Liu>

"""Module description:

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
import numpy as np
from scipy.stats import norm
cdef class HMMR_EM:
    """ Main HMMR EM class.
    
    """
    cdef:
        # define the values that should be returned
        list fragMeans
        list fragStddevs
        
    cdef __init__ ( self, options, genome ):
        # step1: initialize the values of fragMeans and fragStddevs
        # * Set fragment length distribution parameters. 
        #  * Use inputed values to set initial values, if provided. 
        #  * Else use defaults
        #  */

        # cdef:
        # private double[] weights;
        # private double[] mu;
        # private double[] lamda;
        # private double[] data;
        # private double epsilon = 0.0005;
        # private int maxIter=20;
        # private double jump = 1.5;
        

        self.epsilon = 0.0005
        self.maxIter = 20
        self.jump = 1.5

        self.fragMeans = np.array(options.em_means)
        self.fragStddevs = np.array(options.em_stddev)
        
        # double[] fragMeans = new double[4];
        # double[] fragStddevs = new double[4];
        # double[] mode = new double[4];
        # mode[1] = mode[2] = mode[3] = 2;
        # mode[0]=0.5;

        # if (means != null){
        # 	String[] mu = means.split(",");
        # 	for (int i = 0; i < mu.length;i++){
        # 		fragMeans[i] = Double.parseDouble(mu[i]);
        # 	}
        # } else{
        # 	fragMeans[0] = 50.0;
        # 	fragMeans[1] = 200.0; 
        # 	fragMeans[2] = 400.0;
        # 	fragMeans[3] = 600.0;
        # }
        # if (stddevs != null){
        # 	String[] std = stddevs.split(",");
        # 	for(int i = 0;i < std.length;i++){
        # 		fragStddevs[i] = Double.parseDouble(std[i]);
        # 	}
        # } else{
        # 	for (int i = 0;i < fragStddevs.length;i++){
        # 		fragStddevs[i] = 20.0;
        # 	}
        # }

        # step 2, pull large lengths
        #
        # pullLargeLengths puller = new pullLargeLengths(bam, index, minMapQ, genomeStats,java.util.Arrays.copyOf(fragMeans, 4));
        # double[] lengths = puller.getSampledLengths(10);
        # double[] weights = puller.getWeights();
        # puller = null;
        #
        self.pullLargeLengths(options.bamfile, options.min_map_quality, GENOME, options.em_means)
        # the results from pullLargeLengths will be used to initializae HMMR_EM
        
        # step 3, EM trainer
	# //Perform EM training
	# HMMR_EM em = new HMMR_EM(weights,java.util.Arrays.copyOfRange(fragMeans, 1,4),
	# 		java.util.Arrays.copyOfRange(fragStddevs, 1, 4),lengths);
	# em.learn();
        self.learn()
        self.getMeans()
        self.getLambda()
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


    cdef pullLargeLengths(options.bamfile, options.min_map_quality, GENOME, options.em_means):
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
# 			else if(lengths[i] >= cutone && lengths[i] <= cuttwo)
# 				counter2++;
# 			else
# 				counter3++;
# 		}
# 		weights[0] = (double)counter1/(double)lengths.length; 
# 		weights[1] = (double)counter2/(double)lengths.length;
# 		weights[2] = (double)counter3/(double)lengths.length;
		
# 	}

# }
    # java func takes hmmr_em(weights, fragMeans, fragStddev, lengths) all from pullLargeLengths
    # should call hmmr_em during pullLargelengths since these are all used below
    cdef hmmr_em( self, weights, mu, lam, data ):
        self.weights = weights
        self.mu = mu
        self.lam = lam
        self.data = data

    cdef converges( self, value1, value2, epsilon ):
        if (abs(value1 - value2) <= epsilon):
            converged = True
        else:
            converged = False
        return converged
    
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
    cdef get_weighted_density( self, x, mean, lam, weights ):
        normal_dist_density = norm.pdf(x = x, loc = mean, scale = lam)
        return weights * normal_dist_density
        # return np.multiply(weight, normal_dist_density) # make sure both are np.array types

    cdef return_greater( self, data ):
        largest_index = -1
        greatest_value = -1.0
        for i in range(0, len(data)):
            if data[i] > greatest_value:
                greatest_value = data[i]
                largest_index = i
        for i in range(0, len(data)):
            if i != largest_index:
                if data[i] == greatest_value:
                    largest_index = -1
        return largest_index
    

    # https://commons.apache.org/proper/commons-math/javadocs/api-3.1/org/apache/commons/math3/stat/descriptive/moment/Mean.html
    # java version updates mean and std incrementally
    cdef iterate(self):
        temp = []
        counter = []
        total = 0.0
        means = []
        std = []
        for i in range(0, len(hmmr_em.data)):
            for j in range(0. len(hmmr_em.mu)):
                temp[j] = get_weighted_density(hmmr_em.data[i], hmmr_em.mu[j], hmmr_em.lam[j], hmmr_em.weights[j])
            index = return_greater(temp)
            # is this too large of a file to save means,stds for each iteration?
            if index != -1: # if greatest value is not in the last index:
                means[index] = np.mean(hmmr_em.data[0:i]) #check - do we want means[index] or means[i] or means.append()
                stds[index] = np.std(hmmr_em.data[0:i])
                counter[index] += 1.0 # check on this - do we want it to be a list with increasing values?
                total += 1.0
        for i in range(0, len(hmmr_em.mu)): # what are the lengths of means, stds, data, mu?
            hmmr_em.mu[i] = hmmr_em.mu[i] + (jump * (means[i] - hmmr_em.mu[i]))
            hmmr_em.lam[i] = hmmr_em.lam[i] + (jump * (stds[i] - hmmr_em.lam[i]))
            hmmr_em.weights[i] = counter[i] / total

    cdef learn(self):
        converged = False
        itr = 0
        
        while converged == False:
            for a in range(0, len(hmmr_em.mu)):
                old_mu[a] = hmmr_em.mu[a]
                old_lam[a] = hmmr_em.lam[a]
                old_weights[a] = hmmr_em.weights[a]

            iterate()

            counter = 0
            for a in range(0, len(hmmr_em.mu)):
                if converges(old_mu[a], hmmr_em.mu[a], self.epsilon) and converges(old_weights[a], hmmr_em.weights[a], self.epsilon) and converges(old_lam[a], hmmr_em.lam[a], self.epsilon)
                    counter += 1
            if counter == len(hmmr_em.mu):
                converged = True
            itr += 1
            if itr >= self.maxIter:
                break 



            
        # learn function is the major function in HMMR_EM
# 	/**
# 	 * Iterate through the EM algorithm until data convergence
# 	 */
# 	public void learn(){
		
# 		double[] oldMu = new double[mu.length];
# 		double[] oldWeights = new double[mu.length];
# 		double[] oldLam = new double[mu.length];
# 		boolean converged = false;
# 		int iter = 0;
# 		while(!converged){
			
# 			for (int a = 0;a < mu.length;a++){
# 				oldMu[a] = mu[a];
# 				oldLam[a] = lamda[a];
# 				oldWeights[a] = weights[a];
# 			}
			
			
# 			iterate();
			
# 			int counter = 0;
# 			for (int a = 0;a < mu.length;a++){
				
# 				if(converges(oldMu[a],mu[a],epsilon) && converges(oldWeights[a],weights[a],epsilon)
# 						&& converges(oldLam[a],lamda[a],epsilon)){
# 					counter+=1;
# 				}
# 			}
# 			if (counter == mu.length){
				
# 				converged=true;
# 			}
# 			iter+=1;
# 			if (iter >= maxIter){
# 				break;
# 			}
# 			//Output values during iterations
# 			/*
# 			for (int a = 0;a < mu.length;a++){
# 				System.out.println(mu[a]);
				
# 				System.out.println(lamda[a]);
# 				System.out.println(weights[a]);
# 				System.out.println(iter);
# 			}
# 			*/
# 			//System.out.println(iter);
# 		}
		
# 	}
# 	/**
# 	 * Determine if the EM algorithm has converged
# 	 * @param value1 a double representing the value before EM algorithm
# 	 * @param value2 a double representing the value after EM algorithm
# 	 * @param epsilon a double representing the maximum difference allowed to be considered converged
# 	 * @return a boolean indicating whether the values have converged
# 	 */
# 	private boolean converges(double value1,double value2, double epsilon){
# 		if (Math.abs(value1 - value2) <= epsilon){
# 			return true;
# 		}
# 		else{
# 			return false;
# 		}
# 	}
# 	/**
# 	 * Access the weighted density of the multi-variate distribution
# 	 * @param x a double representing the read length
# 	 * @param mean a double representing the distribution mean
# 	 * @param lamda	a double representing the distribution standard deviation
# 	 * @param weight a double representing the distribution weight
# 	 * @return a double representing the weighted density
# 	 */
# 	private double getWeightedDensity(double x, double mean, double lamda,double weight){
		
# 		NormalDistribution dis = new NormalDistribution(mean,lamda);
# 		return weight * dis.density(x);
		
# 	}
# 	/**
# 	 * Return the index of the largest value in an array of doubles
# 	 * @param data an Array of doubles 
# 	 * @return an integer representing thre index of the largest value in the inputted array
# 	 */
# 	private int returnGreater(double[] data){
# 		int largest = -1;
# 		double greatest = -1.0;
# 		for(int i = 0;i < data.length;i++){
# 			if (data[i] > greatest){
# 				greatest = data[i];
# 				largest = i;
# 			}
# 		}
# 		for (int i = 0;i < data.length;i++){
# 			if (i != largest){
# 				if(data[i] == greatest){
# 					largest = -1;
# 				}
# 			}
# 		}
		
# 		return largest;
# 	}
		
# }
