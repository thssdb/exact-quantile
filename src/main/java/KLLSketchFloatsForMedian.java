import org.apache.datasketches.kll.KllFloatsSketch;

import java.util.ArrayList;
import java.util.List;

public class KLLSketchFloatsForMedian {
  KllFloatsSketch sketch;
  final int KllK = 30000;
  //  final double KllEpsilon;
  public KLLSketchFloatsForMedian() {

    sketch = KllFloatsSketch.newHeapInstance(KllK);
//    KllEpsilon = sketch.getNormalizedRankError(false)*1.4;
  }

  public void add(final float value){
    sketch.update(value);
  }

  public List<Long> findResultRange(long K1, long K2) {
    List<Long> result = new ArrayList<>(4);
    long n = sketch.getN();
    result.add((long)sketch.getQuantileLowerBound(1.0d*K1/n));
    result.add((long)sketch.getQuantileUpperBound(1.0d*K2/n));
    return result;
  }
  public boolean isExactResult(){
    return 1.0/sketch.getN()>=sketch.getNormalizedRankError(false)*2;
  }

  public void reset(){
    sketch = KllFloatsSketch.newHeapInstance(KllK);
  }
}