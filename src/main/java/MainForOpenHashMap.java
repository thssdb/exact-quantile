import java.util.Date;
import java.util.List;
import java.util.Random;

public class MainForOpenHashMap {
    int N = 20000000;
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
//        if (i <= KK * 2) num = (i & 1) + (1L << ((i - 1) / 2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/
//        else num = (Math.abs(num) % 65530) << KK;
//        if (i <= KK * 65536)
//            num = (long) ((i - 1) % 65536) << ((i - 1) / 65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/
//        else num = (Math.abs(num) % num_mask) << (KK);
//        num = Math.round(2e18 + 9e18 * random.nextGaussian());
//        if (i % 65500 == 0) random.setSeed(233);
//        num &= 0xFFFFL;
//        num &= 0x000F0F00000F00F0L;
//        num = (num & 65535);

        num = (i & 65535);
        if (i + 65536 * 4 >= N) {
            long tmp = (N - i - 1) / 65536,
                tmp2 = (N - i - 1) % 65536;
            num = tmp2 << (16 * tmp);
        }
        return num;
    }
    public void test0(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        HashMapCacheForQuantile hashMap = new HashMapCacheForQuantile(64,16);
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
            num = nextNumber(num,i);
            hashMap.insert(num,freq);
        }
        List<Long> result = hashMap.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
//        System.out.println(hashMap.DEBUG);
//        System.out.println(hashMap.getRemainingBits());
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test1(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        EclipseCollectionsHashMapForQuantile hashMap = new EclipseCollectionsHashMapForQuantile(64,16);
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
            num = nextNumber(num,i);
            hashMap.insert(num,freq);
        }
        List<Long> result = hashMap.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
//        System.out.println(hashMap.DEBUG);
//        System.out.println(hashMap.getRemainingBits());
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test2(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        HppcHashMapForQuantile hashMap = new HppcHashMapForQuantile(64,16);
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
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test3(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        FastUtilHashMapForQuantile hashMap = new FastUtilHashMapForQuantile(64,16);
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
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test4(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        HashMapOpenCacheForQuantile hashMap = new HashMapOpenCacheForQuantile(64,16);
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
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }

    public static void main(String[] args) {

        System.gc();
        try {
            Thread.sleep(500);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        int TEST_CASE = 3, sum;
        long START_TIME = new Date().getTime();

        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            System.gc();
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            START_TIME = new Date().getTime();
            MainForOpenHashMap main = new MainForOpenHashMap();
            main.test0();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("linked: "+sum/TEST_CASE+"ms");




        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            System.gc();
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            START_TIME = new Date().getTime();
            MainForOpenHashMap main = new MainForOpenHashMap();
            main.test1();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("gc collection: "+sum/TEST_CASE+"ms");




        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            System.gc();
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            START_TIME = new Date().getTime();
            MainForOpenHashMap main = new MainForOpenHashMap();
            main.test2();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("hppc: "+sum/TEST_CASE+"ms");




        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            System.gc();
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            START_TIME = new Date().getTime();
            MainForOpenHashMap main = new MainForOpenHashMap();
            main.test3();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("FastUtil: "+sum/TEST_CASE+"ms");




        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            System.gc();
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            START_TIME = new Date().getTime();
            MainForOpenHashMap main = new MainForOpenHashMap();
            main.test4();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("OpenCache: "+sum/TEST_CASE+"ms");
    }
}
