import java.util.Date;
import java.util.List;
import java.util.Random;

public class MainForKLL {
    static int N = 20000000;
    final long KK = 47, num_mask = 65536;

    long seed;
    Random random = new Random(233);
    private long nextLong(){
        seed ^= (seed << 21);
        seed ^= (seed >>> 35);
        seed ^= (seed << 4);
        return seed;
    }
    private void reset(){seed=2333L;random.setSeed(233);}
    private long nextNumber(long x, int i) {
        long num=x;
        num=x&262143;
        if((i&15)<=3)
            num = nextLong();
        return num;
    }
    public void testAmortized(){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        int com = 6708;
        AmortizedForMedian worker = new AmortizedForMedian(com,com);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=/*Math.abs(nextLong())%23+*/1;
            num = nextNumber(num,i);
            worker.add(num);
        }
        List<Long> result = worker.findResultRange(N/2, N/2+1);
        System.out.println("\t\t\t\t"+result.toString());
//        System.out.println(hashMap.DEBUG);
//        System.out.println(hashMap.getRemainingBits());
//        System.out.println("   end mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void testKLL(){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        KLLSketchDoublesForMedian worker = new KLLSketchDoublesForMedian();
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong();
            num = nextNumber(num,i);
            worker.add(num);
        }
        List<Long> result = worker.findResultRange(N/2, N/2+1);
        System.out.println("\t\t\t\t"+result.toString());
        System.out.println("\t\t"+worker.sketch.getSerializedSizeBytes());
//        System.out.println(hashMap.DEBUG);
//        System.out.println(hashMap.getRemainingBits());
//        System.out.println("   end mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }

    public static void main(String[] args) {

        System.gc();
        try {
            Thread.sleep(500);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        int TEST_CASE = 1, sum;
        long START_TIME;

        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            START_TIME = new Date().getTime();
            MainForKLL main = new MainForKLL();
            main.testKLL();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("KLL Float: "+sum/TEST_CASE+"ms");


        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            START_TIME = new Date().getTime();
            MainForKLL main = new MainForKLL();
            main.testAmortized();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("Amortized: "+sum/TEST_CASE+"ms");



//
//        System.gc();
//        try {
//            Thread.sleep(500);
//        } catch (InterruptedException e) {
//            e.printStackTrace();
//        }
//
//
//
//        TEST_CASE = 40;
//        long start_mem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
//
//        System.out.println("   start mem:"+1.0*start_mem/1024/1024);
////        MainForGCMem[] mainList = new MainForGCMem[TEST_CASE];
//        KLLSketchForMedian[] KLLSketchForMedianList = new KLLSketchForMedian[TEST_CASE];
//        for(int i=0;i<TEST_CASE;i++) {
//            KLLSketchForMedianList[i] = new KLLSketchForMedian();
//            for(int j=0;j<N;j++)KLLSketchForMedianList[i].add(j*1.03f);
////            System.out.println("\t\t\t"+KLLSketchForMedianList[i].sketch.getSerializedSizeBytes());
//        }
////        System.gc();
////        try {
////            Thread.sleep(500);
////        } catch (InterruptedException e) {
////            e.printStackTrace();
////        }
//        long end_mem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
//        System.out.println("   tot mem:"+1.0*(end_mem-start_mem)/1024/1024);
//        System.out.println("   avg mem:"+1.0*(end_mem-start_mem)/TEST_CASE/1024/1024);
//



    }
}
