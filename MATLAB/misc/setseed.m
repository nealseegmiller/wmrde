function setseed(seed)
% set seed for global random stream

stream = RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(stream);