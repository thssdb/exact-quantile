import java.util.Date;
import java.util.List;
import java.util.Random;

public class MainForHashMap {
    int N = 40000000;
    final long KK = 47, num_mask = 65536;

    long seed;
    private long nextLong(){
        seed ^= (seed << 21);
        seed ^= (seed >>> 35);
        seed ^= (seed << 4);
        return seed;
    }
    private void reset(){seed=2333L;}

    public void test1(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        FixedTreap_3 treap_3 = new FixedTreap_3(64,16);
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
//            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<KK;
//            if(i<=KK*65536)num= (long) ((i - 1) % 65536) <<((i-1)/65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<(KK);
//            num&=0x7FFFFFFFFFFFFFFFL;
            treap_3.insert(num,freq);
        }
        List<Long> result = treap_3.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
//        System.out.println(treap_3.DEBUG);
//        System.out.println(treap_3.getRemainingBits());
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test2(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        HashMapForQuantile hashMap = new HashMapForQuantile(64,16);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
//            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<KK;
//            if(i<=KK*65536)num= (long) ((i - 1) % 65536) <<((i-1)/65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<(KK);
//            num = Math.round(2e18 + 9e18 * random.nextGaussian());
//            if(i%65536==0)random.setSeed(233);
//            num&=0x7FFFFFFFFFFFFFFFL;
//            num&=0xFFFF000000000000L;
            hashMap.insert(num,freq);
        }
        List<Long> result = hashMap.findResultIndex(233, 190000000);
        System.out.println(result.toString());
//        System.out.println(hashMap.DEBUG);
        System.out.println(hashMap.getRemainingBits());
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test3(){
        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        HppcHashMapForQuantile hashMap = new HppcHashMapForQuantile(64,16);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
//            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<KK;
//            if(i<=KK*65536)num= (long) ((i - 1) % 65536) <<((i-1)/65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<(KK);
//            num = Math.round(2e18 + 9e18 * random.nextGaussian());
//            if(i%65536==0)random.setSeed(233);
//            num&=0x7FFFFFFFFFFFFFFFL;
            num=(num&65535)%65535;
            hashMap.insert(num,freq);
        }
        List<Long> result = hashMap.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
//        System.out.println(hashMap.DEBUG);
        System.out.println(hashMap.getRemainingBits());
        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test4(){
        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        HashMapCacheForQuantile hashMap = new HashMapCacheForQuantile(64,16);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
//            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<KK;
//            if(i<=KK*65536)num= (long) ((i - 1) % 65536) <<((i-1)/65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<(KK);
//            num = Math.round(2e18 + 9e18 * random.nextGaussian());
//            if(i%65536==0)random.setSeed(233);
//            num&=0x7FFFFFFFFFFFFFFFL;
            num=(num&65535)%65535;
            hashMap.insert(num,freq);
        }
        List<Long> result = hashMap.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
//        System.out.println(hashMap.DEBUG);
        System.out.println(hashMap.getRemainingBits());
        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test5(){
        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        EclipseCollectionsHashMapForQuantile hashMap = new EclipseCollectionsHashMapForQuantile(64,16);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
//            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<KK;
//            if(i<=KK*65536)num= (long) ((i - 1) % 65536) <<((i-1)/65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<(KK);
//            num = Math.round(2e18 + 9e18 * random.nextGaussian());
//            if(i%65536==0)random.setSeed(233);
//            num&=0x7FFFFFFFFFFFFFFFL;
            num=(num&65535)%65535;
            hashMap.insert(num,freq);
        }
        List<Long> result = hashMap.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
//        System.out.println(hashMap.DEBUG);
        System.out.println(hashMap.getRemainingBits());
        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void testMem(){
        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        int N = 3000000;
        Random random = new Random();
        int[] arr = new int[N];
        for(int i=0;i<N;i++)arr[i]=random.nextInt();
        System.out.println("\t   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        long mmp=0;
        for(int i=0;i<N;i++)mmp+=arr[i]^arr[i^1];
        System.out.println(mmp);
//        arr=null;
    }

    public static void main(String[] args) {

        int TEST_CASE = 1, sum;
        long START_TIME = new Date().getTime();

//        System.gc();
//        try {
//            Thread.sleep(500);
//        } catch (InterruptedException e) {
//            e.printStackTrace();
//        }
//        START_TIME = new Date().getTime();
//        for(int i=0;i<TEST_CASE;i++) {
//            MainForHashMap main = new MainForHashMap();
//            main.test1();
//        }
//        System.out.println((new Date().getTime() - START_TIME)/TEST_CASE+"ms");


//        sum=0;
//        for(int i=0;i<TEST_CASE;i++) {
////            System.gc();
////            try {
////                Thread.sleep(500);
////            } catch (InterruptedException e) {
////                e.printStackTrace();
////            }
//            START_TIME = new Date().getTime();
//            MainForHashMap main = new MainForHashMap();
//            main.test2();
//            sum+=(new Date().getTime() - START_TIME);
//        }
//        System.out.println(sum/TEST_CASE+"ms");




        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
//            System.gc();
//            try {
//                Thread.sleep(500);
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
            START_TIME = new Date().getTime();
            MainForHashMap main = new MainForHashMap();
            main.test3();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println(sum/TEST_CASE+"ms");




        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
//            System.gc();
//            try {
//                Thread.sleep(500);
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
            START_TIME = new Date().getTime();
            MainForHashMap main = new MainForHashMap();
            main.test4();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println(sum/TEST_CASE+"ms");




        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
//            System.gc();
//            try {
//                Thread.sleep(500);
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
            START_TIME = new Date().getTime();
            MainForHashMap main = new MainForHashMap();
            main.test5();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println(sum/TEST_CASE+"ms");
    }
}
