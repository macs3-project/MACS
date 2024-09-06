# `hmmratac` HMM model json file format

The HMM trained from `hmmratac` can be saved in a JSON format file and
be loaded later. This option enables users to reuse an ideal hidden
markove model, that captures the signals and relationships among open
chromatin regions, nucleosomes, and backgrounds, from a good quality
data on other ATAC-seq dataset. The JSON data from `hmmratac` is a
JSON dictionary containing:

1. `hmm_type` - either 'gaussian' or 'poisson' for the emission model
2. `hmm_binsize` - the bin size in basepair, used to sample the signals
   across the genome. 
3. `n_features` - this is fixed at 4. In `hmmratac`, the features used
   to train the HMM is the short fragment, the mono-nucleosomal, the
   di-nucleosomal, and the tri-nucleosomal signals.
4. `i_open_region`/`i_nucleosomal_region`/`i_background_region` -
   index number of the three states, including the open region, the
   nucleosomal region and the background region, in the emission and
   transition matrix data, starting from 0.
5. `startprob` - a list of the initial probabilities of the HMM states
   at the first bin of a candidate region for decoding. Ideally, the
   first bin should be more likely a state for the background region,
   and less likely a state for the open region. Check the index
   numbers in `i_open_region`, `i_nucleosomal_region`, or
   `i_background_region` to figure out which value corresponds to
   which state.
6. `transmat` - the transition matrix (list of lists) indicating the
   probabilities that a state can transit to another/or the same
   state. It will be always a 3x3 matrix. If you want to figure the
   transition probability from the state A to the state B, you need to
   find the index number of state A and B, then identify the
   `i_state_A` list, then the `i_state_B` number.
7. `lambda` - this is only available if the model is 'poisson',
   containing the lambda values of Possion models. This represents the
   emission model of each of the three states, so it's a
   3(states)x4(features) matrix.  Check the index numbers in
   `i_open_region`, `i_nucleosomal_region`, or `i_background_region`
   to figure out which list corresponds to which state.
8. `covariance_type` - this is only available while the model type is  
   'gaussian', and is always 'full' currently. 'full' means each state  
   uses a full covariance matrix. 
9. `mean` and `covars`- these are only available if the model is
   'gaussian'. If the emission follows Gaussian model, each state will
   have a mean value and a 4x4 matrix for the full covariance matrix.
   Check the index numbers in `i_open_region`, `i_nucleosomal_region`,
   or `i_background_region` to figure out which mean value or 4x4
   matrix correspond to which state.
