import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.apache.datasketches.kll.KllDoublesSketch;

public class DataSketchesKLLForExact {
  public KllDoublesSketch sketch;
  int KLLK = 0;
  //  final double KllEpsilon;
  public DataSketchesKLLForExact(int KLL_K) {
    KLLK = KLL_K;
    sketch = KllDoublesSketch.newHeapInstance(KLLK);
//    KllEpsilon = sketch.getNormalizedRankError(false)*1.4;
  }

  public void update(final double value){
    sketch.update(value);
  }

  public double[] findResultRange(long K1, long K2) {
    DoubleArrayList result = new DoubleArrayList(4);
    long n = sketch.getN();
    if(sketch.getN()<=KLLK){
      result.add(sketch.getQuantile(1.0d * K1 / n));
      result.add(sketch.getQuantile(1.0d * K2 / n));
      result.add(-233);
//      System.out.println("\t!!!OVER EXACT GOT");
    }else {
      result.add(sketch.getQuantileLowerBound(1.0d * K1 / n));
      result.add(sketch.getQuantileUpperBound(1.0d * K2 / n));
    }
//    // rank of x is the sum of weight of items < x in datasketches kll.
//    long nextK1 = K1-Math.round(sketch.getRank(result.getDouble(0))*n),nextK2 = K2-(K1-nextK1);
//    result.add(nextK1);
//    result.add(nextK2);
    return result.toDoubleArray();
  }

  public void reset(){
    sketch = KllDoublesSketch.newHeapInstance(KLLK);
  }

  static int calcKUnderLimitedMemory(int maxMemoryByte,long queryN){
    int maxItems = maxMemoryByte/8;
    int K=10;
    for(int addK=Integer.highestOneBit(maxItems);addK>=1;addK>>=1){
      int cntK=K+addK;
      KllDoublesSketch tmpSketch=KllDoublesSketch.newHeapInstance(cntK);
      int cntBytes=0;
      for(int i=0;i<queryN;i++){
        tmpSketch.update(i);
      }
      if(tmpSketch.getCurrentUpdatableSerializedSizeBytes()<=maxMemoryByte)
        K+=addK;
    }
    return K;
  }

}