/*
 * =====================================================================================
 *
 *       Filename:  Constants.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2014年11月02日 00时33分07秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xu Fei (Fly), xufei@lsec.cc.ac.cn
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <malloc.h>
#include "Constants.h"
//对向量a的 a[left:right]进行排练, 数组长度 right-left+1

void QuickSort_Int(int a[],int left,int right)
{
  int i,j,temp;
  i=left;
  j=right;
  //temp=a[left];
  if(left>right)
    return;
  temp=a[left];
  while(i!=j)
   {
     while(a[j]>=temp&&j>i)
       j--;
     if(j>i)
       a[i++]=a[j];
     while(a[i]<=temp&&j>i)
       i++;
     if(j>i)
       a[j--]=a[i];     
  }
  a[i]=temp;
  QuickSort_Int(a,left,i-1);  
  QuickSort_Int(a,i+1,right);
}


//对向量a的 a[left:right]进行排练, 数组长度 right-left+1
void QuickSort_Double(double a[],int left,int right)
{
  int i,j;
  double temp;
  i=left;
  j=right;
  if(left>right)
    return;
  temp=a[left];
  while(i!=j)
   {
     while(a[j]>=temp&&j>i)
       j--;
     if(j>i)
       a[i++]=a[j];
     while(a[i]<=temp&&j>i)
       i++;
     if(j>i)
       a[j--]=a[i];     
  }
  a[i]=temp;
  QuickSort_Double(a,left,i-1);  
  QuickSort_Double(a,i+1,right);
}


double GetTime()
{
  struct rusage usage;
  double ret;
    
  if(getrusage(RUSAGE_SELF, &usage) == -1) 
    printf("Error in GetTime!\n");

  ret = ((double) usage.ru_utime.tv_usec)/1000000;

  ret += usage.ru_utime.tv_sec;

  return ret;
}


double GetMemory()
{

  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();

  return (double)(MALLINFO.usmblks+MALLINFO.uordblks)/1048576;  

}

