/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

import java.util.ArrayList;
import java.util.List;

/**
 * This is a bitwise trie (01 trie) with a maximum size.
 */
public class FixedTreap {
  private static final int maxSize = (1 << 16) + 1;
  //  private static final int SHR = 1;
  private long deltaForUnsignedCompare;
  int[] tmpQueue;
  int tmpQueueLength;
  int[] tmpStack;
  int tmpStackTop;
  int[] trashStack;
  int trashTop;
  int[] lc, rc, fa, rnd;
  long[] value, count;
  int root;
  long randSeed;
  int bitsOfValue, remainingBits;
  long DEBUG=0;


  public FixedTreap(int bits) { // TODO use bucket when remainingBits is too small(16)
    tmpQueue = new int[maxSize + 1];
    tmpStack = new int[233];
    trashStack = new int[maxSize + 1];
    trashTop = maxSize;
    for (int i = 1; i <= maxSize; i++)
      trashStack[i] = maxSize-i+1;
    lc = new int[maxSize + 1];
    rc = new int[maxSize + 1];
    fa = new int[maxSize + 1];
    rnd = new int[maxSize + 1];
    value = new long[maxSize + 1];
    count = new long[maxSize + 1];
    root = 0;
    randSeed = System.currentTimeMillis();
    bitsOfValue = remainingBits = bits;
    if (bits == 64) deltaForUnsignedCompare = 1L << 63; // unsigned long
    else deltaForUnsignedCompare = 0;
  }

  private int getNewNode() {
    int x = trashStack[trashTop--];
    lc[x] = rc[x] = fa[x] = 0;
    count[x] = 0;
    return x;
  }

  private void freeNode(int x) {
    trashStack[++trashTop] = x;
  }
  private int nextRandInt(){ // xorShift
    randSeed ^= (randSeed << 21);
    randSeed ^= (randSeed >>> 35);
    randSeed ^= (randSeed << 4);
    return (int)randSeed;
  }

  private int insert(int x, long num, long freq) {
    if (x == 0) {
      x = getNewNode();
      rnd[x] = nextRandInt();
      value[x] = num;
      count[x] = freq;
    } else if (num == value[x]) {
      count[x] += freq;
    } else if ((num ^ deltaForUnsignedCompare) < (value[x] ^ deltaForUnsignedCompare)) {
      int l = lc[x] = insert(lc[x], num, freq);
      fa[l] = x;
      if (rnd[l] < rnd[x]) { // rotate
        lc[x] = rc[l];
        fa[rc[l]] = x;
        rc[l] = x;
        fa[l] = fa[x];
        fa[x] = l;
        return l;
      }
    } else if ((num ^ deltaForUnsignedCompare) > (value[x] ^ deltaForUnsignedCompare)) {
      int r = rc[x] = insert(rc[x], num, freq);
      fa[r] = x;
      if (rnd[r] < rnd[x]) { // rotate
        rc[x] = lc[r];
        fa[lc[r]] = x;
        lc[r] = x;
        fa[r] = fa[x];
        fa[x] = r;
        return r;
      }
    }
    return x;
  }

  private void deleteNode01(int x) { // x won't be root; x has less than 2 child
    int ch = lc[x] + rc[x];
    fa[ch] = fa[x];
//    if(lc[fa[x]]!=x&&rc[fa[x]]!=x)System.out.println("[DEBUG]????");
    if (lc[fa[x]] == x) lc[fa[x]] = ch;
    else rc[fa[x]] = ch;
    freeNode(x);
  }
//
//  private int findSHR() {
//    int suffixLen = 62;
//    long suffixMask = ((1L << suffixLen) - 1);
//    long suffix = value[tmpQueue[1]] & suffixMask;
//    long val;
//    for (int i = 2; suffixLen > 0 && i <= tmpQueueLength; i++) {
//      val = value[tmpQueue[i]];
//      while (suffixLen > 0 && (val & suffixMask) != suffix) {
//        suffixLen -= 1;
//        suffixMask >>>= 1;
//        suffix &= suffixMask;
//      }
//    }
//    System.out.println("[DEBUG treap] findSHR  suffixLEN:"+suffixLen+"   suffix:"+suffix);
//    return suffixLen + 1;
//  }

  private int findSHR() {
    long x,y,xor;
    int SHR=remainingBits-16, tmpLen;
    for (int i = 2; i <= tmpQueueLength; i++) {
      x = value[tmpQueue[i-1]];
      y = value[tmpQueue[i]];
      xor = x^y;
      if((xor>>>SHR)>0)continue;
      tmpLen = 0;
//      System.out.println("[DEBUG treap] findSHR  xor:"+xor);
      for(int j = 32;j>0;j>>=1)
        if((xor>>>j)>0){
          tmpLen+=j;
          xor>>>=j;
        }
      tmpLen+=xor;
      SHR = Math.min(SHR, tmpLen);
    }
//    System.out.println("[DEBUG treap] findSHR  :"+SHR);
    return SHR;
  }
  // TODO insert非递归化；rebuild重构(rnd swap)
  private void rebuild() {
    tmpQueueLength = 0;
    dfs(root);
    int SHR = findSHR();
    remainingBits -= SHR;
    deltaForUnsignedCompare = 0;
    for (int i = 1; i <= tmpQueueLength; i++)
      value[tmpQueue[i]] >>>= SHR;

    int x, preX = tmpQueue[1];
    for (int i = 2; i <= tmpQueueLength; i++) {
      x = tmpQueue[i];
      if (value[x] == value[preX]) {
        if (lc[x] != 0 && rc[preX] == 0) {
          deleteNode01(preX);
          count[x] += count[preX];
          preX = x;
        } else { // lc[x]==0&&rc[preX]!=0
          deleteNode01(x);
          count[preX] += count[x];
        }
      } else preX = x;
    }
  }

  private void dfs(int x) {
//    System.out.println("[DEBUG TREAP] dfs  x:"+x+"  lc:"+lc[x]+"  rc:"+rc[x]+"     "+value[x]+"  "+count[x]);
    if (lc[x] != 0) dfs(lc[x]);
    tmpQueue[++tmpQueueLength] = x;
    if (rc[x] != 0) dfs(rc[x]);
  }
  private void dfs_DEBUG(int x) {
    System.out.println("[DEBUG TREAP] dfs  x:"+x+"  lc:"+lc[x]+"  rc:"+rc[x]+"     "+value[x]+"  "+count[x]+"\t\t"+rnd[x]);
    if (lc[x] != 0) dfs_DEBUG(lc[x]);
    if (rc[x] != 0) dfs_DEBUG(rc[x]);
  }

  private int getSize() {
    return maxSize - trashTop;
  }


  public void insert(long num, long freq) {
//    System.out.println("[DEBUG TREAP] insert  num:"+num+"  freq:"+freq);
    num>>>=bitsOfValue-remainingBits;
//    System.out.println("[DEBUG TREAP] insert  num:"+num+"  freq:"+freq);
    root = insert(root, num, freq);
    while (trashTop == 0) { // when size >= 2^16+1
      rebuild();
    }
//    dfs_DEBUG(root);System.out.println("\n\n");
  }

  public List<Long> findResultIndex(long K1, long K2) {
    List<Long> result = new ArrayList<>(8);
    long sum = 0;
    tmpQueueLength = 0;
    dfs(root);
    for (int i = 1; i <= tmpQueueLength; i++) {
//      System.out.println(count[tmpQueue[i]] + "  " + value[tmpQueue[i]]);
      sum += count[tmpQueue[i]];
      if (sum >= K1 && result.size() == 0) {
        result.add(value[tmpQueue[i]]);
        result.add(sum - count[tmpQueue[i]]);
      }
      if (sum >= K2 && result.size() == 2) {
        result.add(value[tmpQueue[i]]);
        result.add(sum - count[tmpQueue[i]]);
        break;
      }
    }
    return result;
  }

  public int getRemainingBits() {
    return remainingBits;
  }
}
