import org.apache.datasketches.kll.KllDoublesSketch;

import java.util.Arrays;
import java.util.Random;

public class MainForKLLStatistics {
    static int N = 20000000;
    final long KK = 47, num_mask = 65536;

    long seed;
    Random random = new Random(233);
    private long nextLong(){
        return random.nextLong();
//        seed ^= (seed << 21);
//        seed ^= (seed >>> 35);
//        seed ^= (seed << 4);
//        return seed;
    }
    private void reset(){seed=2333L;random.setSeed(233);}
    private long nextNumber(long x, int i) {
        long num=x;
        num=x&262143;
        if((i&15)<=3)
            num = nextLong();
        return num;
    }

    public void testKLL_PAGE(){
        final int NNN_PAGE=7000;
        System.out.println("\t\t KKK_K8:"+8+"  size:"+
            KllDoublesSketch.getMaxSerializedSizeBytes(8, NNN_PAGE, false));
        KllDoublesSketch pageSketch = KllDoublesSketch.newHeapInstance(NNN_PAGE);
        int KKK_PAGE=0;

        for(;KKK_PAGE<=2000;KKK_PAGE++) {
            int size = KllDoublesSketch.getMaxSerializedSizeBytes(KKK_PAGE, NNN_PAGE, false);
            if (size>=1024)break;
        }KKK_PAGE-=1;
        System.out.println("\t\t KKK_001:"+KKK_PAGE+"  size:"+
            KllDoublesSketch.getMaxSerializedSizeBytes(KKK_PAGE, NNN_PAGE, false));
        KllDoublesSketch sketch_001 = KllDoublesSketch.newHeapInstance(KKK_PAGE);
        System.out.println(sketch_001.getNormalizedRankError(false));

        for(KKK_PAGE=0;KKK_PAGE<=2000;KKK_PAGE++) {
            int size =KllDoublesSketch.getMaxSerializedSizeBytes(KKK_PAGE, NNN_PAGE, false);
            if (size>=10240)break;
        }KKK_PAGE-=1;
        System.out.println("\t\t KKK_01:"+KKK_PAGE+"  size:"+
            KllDoublesSketch.getMaxSerializedSizeBytes(KKK_PAGE, NNN_PAGE, false));
        KllDoublesSketch sketch_01 = KllDoublesSketch.newHeapInstance(KKK_PAGE);
        System.out.println(sketch_01.getNormalizedRankError(false));
    }
    public void testKLL_CHUNK(){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
//        final int KLL_NUM = 2, KKK = 10700, NNN_CHUNK=1600000;
//        KLLStatistics[] KLLList = new KLLStatistics[KLL_NUM];
//        reset();
//        for(int T=0;T<KLL_NUM;T++) {
//            KLLStatistics worker = new KLLStatistics(KKK);
//            KLLList[T] = worker;
////            System.out.println(KllDoublesSketch.getMaxSerializedSizeBytes(40,7000,false));
//            for(int i=0;i<NNN;i++)
//                worker.add(i);
////            System.out.println(worker.sketch.toString(false,false));
//            System.out.println(worker.sketch.getSerializedSizeBytes());
//            System.out.println(worker.sketch.getNormalizedRankError(false));
//        }
//        System.out.println("MERGE..........");
//        KllDoublesSketch UnionSketch = KllDoublesSketch.newHeapInstance(KKK);
//        for(int i=0;i<KLL_NUM;i++)
//            UnionSketch.merge(KLLList[i].sketch);
//        System.out.println(UnionSketch.getSerializedSizeBytes());
//        System.out.println(UnionSketch.getNormalizedRankError(false));
//        System.out.println(UnionSketch.toString(true,false));
//
//        System.out.println("THEORY..........");
//        System.out.println(KllDoublesSketch.getMaxSerializedSizeBytes(KKK, NNN*KLL_NUM, false));
        final int NNN_CHUNK = 1600000;
        int KKK_CHUNK=0;
        for(KKK_CHUNK=0;KKK_CHUNK<=20000;KKK_CHUNK++) {
            int size =KllDoublesSketch.getMaxSerializedSizeBytes(KKK_CHUNK, NNN_CHUNK, false);
            if (size>=256*1024)break;
        }KKK_CHUNK-=1;
        System.out.println("\t\t KKK_CHUNK:"+KKK_CHUNK+"  size:"+
            KllDoublesSketch.getMaxSerializedSizeBytes(KKK_CHUNK, NNN_CHUNK, false));
        KllDoublesSketch sketch_CHUNK = KllDoublesSketch.newHeapInstance(KKK_CHUNK);
        System.out.println(sketch_CHUNK.getNormalizedRankError(false));
    }
    public void testKLL_ENOUGHK(){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        final int KLL_NUM = 1, KKK = 45000, NNN=40000;
        reset();
        KllDoublesSketch sketch = KllDoublesSketch.newHeapInstance(KKK);
        for(int i=0;i<NNN;i++)sketch.update(2.0*i);
        System.out.println(sketch.getSerializedSizeBytes());
        System.out.println(sketch.getNormalizedRankError(false));
//        System.out.println(sketch.toString(false,false));
        int K1 = (NNN+1)/2-1, K2 = NNN/2;
        System.out.println(sketch.getQuantileLowerBound(1.0*K1/NNN));
        System.out.println(sketch.getQuantileUpperBound(1.0*K1/NNN));
        System.out.println(sketch.getQuantileLowerBound(1.0*K2/NNN));
        System.out.println(sketch.getQuantileUpperBound(1.0*K2/NNN));
        System.out.println(1.0/sketch.getN()>=sketch.getNormalizedRankError(false)*2);
        System.out.println(sketch.toString(true,false));
        System.out.println(sketch.getQuantile(1.0*K1/NNN)+", "+sketch.getQuantile(1.0*K2/NNN));
        System.out.println(Arrays.toString(sketch.getQuantiles(5)));
        System.out.println("N:"+sketch.getN()+"  K:"+sketch.getK());
    }
    public void testKLL_CHUNK_MERGE_PAGE(){
        final int NNN_PAGE=7000,NNN_CHUNK = 1600000,PAGES=10;
        KllDoublesSketch chunkSketch = KllDoublesSketch.newHeapInstance(11363);
        for(int i=1;i<=NNN_CHUNK;i++)chunkSketch.update(1.0*i);
        KllDoublesSketch UnionSketch = KllDoublesSketch.newHeapInstance(11363);
        UnionSketch.merge(chunkSketch);
        System.out.println("\t\t eps1:"+UnionSketch.getNormalizedRankError(false));

        for(int T=0;T<PAGES;T++) {
            KllDoublesSketch pageSketch = KllDoublesSketch.newHeapInstance(32);
            for (int i = 1; i <= NNN_PAGE; i++) pageSketch.update(1.0 * i + T*NNN_PAGE+NNN_CHUNK);
            UnionSketch.merge(pageSketch);
        }
        System.out.println("\t\t eps_merged:"+UnionSketch.getNormalizedRankError(false));
        int NNN = NNN_CHUNK+PAGES*NNN_PAGE;
        int K1 = (NNN+1)/2-1, K2 = NNN/2;
        System.out.println("\t\t v1,v2:"+UnionSketch.getQuantile(1.0*K1/NNN)+","+UnionSketch.getQuantile(1.0*K2/NNN));
        System.out.println("\t\t N,mid:"+NNN+","+NNN/2);
    }

    public static void main(String[] args) {
        MainForKLLStatistics main;
        main = new MainForKLLStatistics();
        main.testKLL_PAGE();System.out.println("\n\n");
//        main.testKLL_CHUNK();System.out.println("\n\n");
//        main.testKLL_CHUNK_MERGE_PAGE();System.out.println("\n\n");
//        main.testKLL_ENOUGHK();
    }
}
