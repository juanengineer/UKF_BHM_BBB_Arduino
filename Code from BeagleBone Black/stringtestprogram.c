#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int main(){

char str1[30] = "Total";
char str2[30] = "Total Battery Voltage = 7";
char str3[25];
char str4[25];

int len; 
int res;

int strpos = 0;
float val;

//strcpy(str3, str1);
//printf("strcpy( str3, str1 ) : %s\n", str3 );

//strcat (str1, str2);
//printf("strcat( str1, str2):     %s\n", str1 );

//len = strlen(str1);
//printf("strlen(str1) : %d\n", len );

res = strcmp(str1, str2);

printf("String comparision: %d\n",res);

if (res == 0) {
str4[0] = str1[24];
printf("The extracted string is: %s\n",str4);

}
else{
printf("The correct string was not identified\n");

}

char *str5 = "rl";
char *ans;
ans = strchr(str5,'l');
printf("The location of str 5 is %p\n", str5);
printf("The location of l is %p\n", ans);

char *str6 = "World";
char *ans1;
/* ans1 = strstr(str6,str5);
ans1 = strstr(str6,"rl");
printf("The location of str6 is %p\n",str6);
printf("The location of rl within str 6 is %p\n",ans1);


strpos = ans1 - str6;
printf("strpos = %d\n",strpos);
*/

ans1 = strstr(str2,"Voltage");
strpos = ans1 - str2;
printf("strpos = %d\n",strpos);
str4[0] = str2[strpos+10]; 
printf("Value is= %s\n",str4);
val = atof(str4);
printf("Convert string to float and then add 2: %f\n",val+2);
return 0;
}
