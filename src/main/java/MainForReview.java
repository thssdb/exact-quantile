import java.util.HashSet;
import java.util.Iterator;

public class MainForReview {
    public static class pp{
        public int s,d,l;
        public pp(int _s,int _d,int _l){s=_s;d=_d;l=_l;}
        public String toString(){
            return "("+s+","+d+","+l+") ";
        }

    }

    public static void main(String[] args) {
        int sigma = 3;
        int n = 600;
        int[] a = new int[n];
        HashSet<pp> P = new HashSet<>();
        HashSet<pp> MP = new HashSet<>();
        for(int i=0;i<n/2;i++)a[i]=i*i;
        for(int i=n/2;i<n-1;i++)a[i]=(n/2-1)*(n/2-1)+i-(n/2-1);
        a[n-1] = (n-1)*(n-1);
        int MAX = a[n-1];
        for(int i=0;i<n;i++){
            int t = a[i];
            HashSet<Integer> ed = new HashSet<>();
            Iterator<pp> itp = P.iterator();
            while(itp.hasNext()) {
                pp p = itp.next();
                if (t - p.s == p.l * p.d) {
                    p.l++;
                    ed.add(p.d);
                } else if (t - p.s > p.l * p.d) {
                    if (p.l >= sigma)
                        MP.add(p);
                    itp.remove();
                }
            }
            for(int j=0;j<i;j++) {
                int s = a[j];
                if ((t - s)*(sigma-2)<=MAX-t && !ed.contains(t-s)){
                    P.add(new pp(s,t-s,2));
                }
            }
//            if(i>=n/2)
//            System.out.println("\t i:"+(i+1) + " + P.size:"+P.size() + "\t\t"+(i*i/3));
//            System.out.println("\t P:"+P);
        }
        for(pp p:P)
            if(p.l>=sigma)
                MP.add(p);
        System.out.println("\t\t\t MP:"+MP);
    }
}
