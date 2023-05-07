import java.util.Date;
import java.util.List;
import java.util.Random;

public class MainForTDigestCached {
    int N = 40000000;
    final long KK = 47, num_mask = 65536;
    final static int com = 6708;

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
//        num^=1L<<63;
//        num = (x&127)+128;
//        if(Math.abs(nextLong())%10000<=4500)
//            num = nextLong();
        return num;
    }
    public void test1(){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
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
    public void test2(){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        TDigestCachedForMedian worker = new TDigestCachedForMedian(com,com);
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

    public static void main(String[] args) {

        System.gc();
        try {
            Thread.sleep(500);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        int TEST_CASE = 5, sum;
        long START_TIME = new Date().getTime();

        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
//            System.gc();
//            try {
//                Thread.sleep(500);
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
            START_TIME = new Date().getTime();
            MainForTDigestCached main = new MainForTDigestCached();
            main.test1();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("amortized: "+sum/TEST_CASE+"ms");

        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
//            System.gc();
//            try {
//                Thread.sleep(500);
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
            START_TIME = new Date().getTime();
            MainForTDigestCached main = new MainForTDigestCached();
            main.test2();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("amortized cached: "+sum/TEST_CASE+"ms");


//
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
//
//        TEST_CASE = 40;
//        long start_mem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
//
//        System.out.println("   start mem:"+1.0*start_mem/1024/1024);
////        MainForGCMem[] mainList = new MainForGCMem[TEST_CASE];
//        AmortizedForMedian[] AmortizedForMedianList = new AmortizedForMedian[TEST_CASE];
//        for(int i=0;i<TEST_CASE;i++) {
//
////            mainList[i] = new MainForGCMem();
////        System.out.println("\t\t   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
//            AmortizedForMedianList[i] = new AmortizedForMedian(com, com);
////        System.out.println("\t\t   end mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
////            mainList[i].test4(GCHashMapForMemList[i]);
//        }
//
//        System.gc();
//        try {
//            Thread.sleep(500);
//        } catch (InterruptedException e) {
//            e.printStackTrace();
//        }
//        long end_mem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
//        System.out.println("   tot mem:"+1.0*(end_mem-start_mem)/1024/1024);
//        System.out.println("   avg mem:"+1.0*(end_mem-start_mem)/TEST_CASE/1024/1024);
//



    }
}
