import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.apache.commons.math3.distribution.LogNormalDistribution;

import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.TreeMap;

public class Main {
    int N = 2000000;
    final long KK = 47, num_mask = 65536;

    long seed;
    private long nextLong(){
        seed ^= (seed << 21);
        seed ^= (seed >>> 35);
        seed ^= (seed << 4);
        return seed;
    }
    private void reset(){seed=2333L;}


    public void test0(){
        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        HashMap<Long,Long> hashMap = new HashMap<>(1<<16, 1.0f);
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=nextLong()%23+1;
            hashMap.merge(num&0xFFFFL, freq, Long::sum);
        }
        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test0_A(){
        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        TreeMap<Long,Long> TreeMap = new TreeMap<>();
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=nextLong()%23+1;
            TreeMap.merge(num&0xFFFFL, freq, Long::sum);
        }
        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }

    public void test1(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        FixedTreap treap = new FixedTreap(64);
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=nextLong()%23+1;
            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));else if(i<=KK*2+65536)num= (long) i <<KK;else num=(Math.abs(num)%num_mask)<<KK;
            treap.insert(num,freq);
//            trie.insert(random.nextLong(), random.nextLong(1, 233));
        }
        List<Long> result = treap.findResultIndex(233, 20000);
        System.out.println(result.toString());
        System.out.println(treap.getRemainingBits());
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test2(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        FixedTreap_2 treap_2 = new FixedTreap_2(64);
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=nextLong()%23+1;
//            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<KK;
            if(i<=KK*65536)num= (long) ((i - 1) % 65536) <<((i-1)/65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<(KK);
            treap_2.insert(num,freq);
//            treap.insert(random.nextLong(), random.nextLong(1, 233));
        }
        List<Long> result = treap_2.findResultIndex(233, 20000);
        System.out.println(result.toString());
//        System.out.println(treap_2.DEBUG);
        System.out.println(treap_2.getRemainingBits());
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test3(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        FixedTreap_3 treap_3 = new FixedTreap_3(64,16);
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
//            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<KK;
            if(i<=KK*65536)num= (long) ((i - 1) % 65536) <<((i-1)/65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<(KK);
//            num&=0x7FFFFFFFFFFFFFFFL;
            treap_3.insert(num,freq);
//            treap.insert(random.nextLong(), random.nextLong(1, 233));
        }
        List<Long> result = treap_3.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
//        System.out.println(treap_3.DEBUG);
//        System.out.println(treap_3.getRemainingBits());
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test4(){
//        System.out.println("   start mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        HashMapForQuantile hashMap = new HashMapForQuantile(64,16);
        reset();
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=Math.abs(nextLong())%23+1;
//            if(i<=KK*2)num=(i&1)+(1L<<((i-1)/2));/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<KK;
            if(i<=KK*65536)num= (long) ((i - 1) % 65536) <<((i-1)/65536);/*else if(i<=KK*2+65536)num= (long) i <<KK;*/else num=(Math.abs(num)%num_mask)<<(KK);
//            num&=0x7FFFFFFFFFFFFFFFL;
            hashMap.insert(num,freq);
//            treap.insert(random.nextLong(), random.nextLong(1, 233));
        }
        List<Long> result = hashMap.findResultIndex(233, 190000000);
//        System.out.println(result.toString());
//        System.out.println(hashMap.DEBUG);
//        System.out.println(hashMap.getRemainingBits());
//        System.out.println("   end mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
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

//    public static void main(String[] args) {
//        int n=(int)2e9;
//        ZipfDistribution zipf = new ZipfDistribution(n/10,0.5);
//        IntArrayList data=new IntArrayList();
//        for(int i=0;i<n;i++)data.add(zipf.sample());
//        data.unstableSort(null);
//        int CNT=0;
//        for (int i = 0; i < n; i++) if(i==0||data.getInt(i)!=data.getInt(i-1))CNT++;
//        int CNT2=0;
////        for (int i = 0; i < n; i++) if(data.getInt(i)<=n/10/2)CNT2++;
//        System.out.println("\t\tCOUNT:\t"+CNT+"\t\t50%_COUNT:\t"+CNT2);
//        System.out.println("\t\tith:\t\t\t0.25,mid,0.75:"+data.getInt(n/4)+","+data.getInt(n/2)+","+data.getInt(n/4*3));
//    }
    //public static void main(String[] args) {
    //    int n=(int)2e9;
    //    UniformRealDistribution R = new UniformRealDistribution(0,1);
    //    Random bit = new Random();
    //    double Emax=300;
    //    DoubleArrayList data = new DoubleArrayList();
    //    for(int i=0;i<(int)1e3;i++){
    //        double x = Math.pow(10, Emax*(Math.pow(R.sample(),2)*2-1));
    //        if(bit.nextBoolean())
    //            x=-x;
    //        data.add(x);
    //    }
    //    data.sort(Double::compare);
    //    System.out.println("\t\t"+data);
    //}
    public static void main(String[] args) {
        int n=(int)1e5;
        LogNormalDistribution R = new LogNormalDistribution(4,0.0004);
        Random bit = new Random();
        double Emax=300;
        DoubleArrayList data = new DoubleArrayList();
        for(int i=0;i<n;i++){
            double x = R.sample();
            data.add(x);
        }
        data.sort(Double::compare);
        System.out.println("min:\t\t"+data.subList(0,10));
        System.out.println("max:\t\t"+data.subList(n-10,n));
        System.out.println("avg:\t\t"+data.getDouble(n/2));
        for(int j=10,i=0;i<j;i++)
            System.out.println("i:"+i+"\t\t"+data.getDouble(n/(j+1)*(i+1)));
        int disC=0;
        for(int i=1;i<n;i++)if(data.getDouble(i-1)!=data.getDouble(i))disC++;
        System.out.println("\t\tdistinct count:\t"+disC);
    }
}
