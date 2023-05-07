import java.util.Date;

public class MainForSimulateKLLCompactNum {
    static int TEST_CASE=100,maxMemoryByte=4096,maxN=25000000;//8192*6771;
    public static void main(String[] args) {

        long START_T = new Date().getTime();
        for(int T=0;T<TEST_CASE;T++){
            KLLSketchLazyEmptyForSimuCompact worker = new KLLSketchLazyEmptyForSimuCompact(maxMemoryByte/8);
//            for(int i=1;i<=10;i++) {
//                int[] compactNum = worker.simulateCompactNumGivenN(maxN*i/10);
////                System.out.println("\tcompactionNum:\t" + Arrays.toString(compactNum));
//            }
            worker.simulateCompactNumGivenN(maxN);
        }
        long ALL_T=(new Date().getTime()-START_T);
        System.out.println("\t\t\tALL_TIME:"+ALL_T+"\t\tavg_t:"+1.0*ALL_T/TEST_CASE);
    }
}
