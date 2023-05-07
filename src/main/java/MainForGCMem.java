import java.util.Date;
import java.util.List;
import java.util.Random;

public class MainForGCMem {
    int N = 2000000;
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
//        num = (num & 65535);
        return num;
    }
    public void test4(GCHashMapForMem hashMap){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
            num = nextNumber(num,i);
            hashMap.insert(num,freq);
        }
        List<Long> result = hashMap.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
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

        long START_TIME = new Date().getTime();
        int TEST_CASE = 40, sum;
        long start_mem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();

        System.out.println("   start mem:"+1.0*start_mem/1024/1024);
//        MainForGCMem[] mainList = new MainForGCMem[TEST_CASE];
        GCHashMapForMem[] GCHashMapForMemList = new GCHashMapForMem[TEST_CASE];
        sum=0;
        for(int i=0;i<TEST_CASE;i++) {

//            mainList[i] = new MainForGCMem();
//        System.out.println("\t\t   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
            GCHashMapForMemList[i] = new GCHashMapForMem(64,16, 1<<16);
//        System.out.println("\t\t   end mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
//            mainList[i].test4(GCHashMapForMemList[i]);
        }

        System.gc();
        try {
            Thread.sleep(2000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        System.out.println("time GCHashMapForMem: "+(new Date().getTime() - START_TIME)/TEST_CASE+"ms");
        long end_mem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        System.out.println("   end mem:"+1.0*(end_mem-start_mem)/1024/1024);
        System.out.println("   avg mem:"+1.0*(end_mem-start_mem)/TEST_CASE/1024/1024);




    }
}
