// inspired by t-Digest by Ted Dunning. See https://github.com/tdunning/t-digest
// This is a simple implementation with radix sort and K0.
// Clusters are NOT strictly in order.

import com.tdunning.math.stats.MergingDigest;
import com.tdunning.math.stats.TDigest;

// 没用。
public class TDigestTedForStatMerge {
  public final long deltaForUnsigned;
  public final int compression;
  public TDigest cluster;
  public int clusterNum, clusterNumMemLimit, clusterNumSeriLimit;
  public int maxSeriByte, maxMemByte;
  public long totN;
  boolean sorted = true;

  public TDigestTedForStatMerge(int maxMemByte, int maxSeriByte) {
    this.maxMemByte = maxMemByte;
    this.maxSeriByte = maxSeriByte;
    this.clusterNumMemLimit = maxMemByte / (2*8);
    clusterNumSeriLimit = maxSeriByte / (2*8);
    clusterNum = 0;
    totN = 0;
    int rate = 1;
    int size = (int)(clusterNumMemLimit/(1+rate*20.0/16));
    int bufferSize = size*rate;
    this.compression = 2*size; //  cluster:buffer = 1:5
    cluster = new MergingDigest(this.compression,bufferSize,size);
    deltaForUnsigned = 1L << (64 - 1);
  }
  public void update(double value) {
    totN++;
    cluster.add(value);
  }

  public double quantile(double q){
    return cluster.quantile(q);
  }




  public void reset() {
    clusterNum = 0;
    totN = 0;
  }

}