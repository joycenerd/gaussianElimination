#  高斯消去法

> 高斯消去法介紹

- 方程組一定是Ax=b，而A是方程組係數所構成的矩陣，b是每個方程式等號右邊的那個常數
- 如果有x1,x2...xn，n個未知數，那A會是n*n的矩陣，而b會是 n*1 的矩陣
- 目標是把A和b合併成擴增矩陣，然後做高斯消去法，把主對角線消為1，其他地方變成0，原本b的那一個column的數字就是未知數的解
- 在算未知數解的時候，在過程中一邊算行列式值
- 有可能無解，就是某一個row全為0，也就是變成有效的方程式小於n也就是未知數的數量，那此方程組無唯一解，那此時行列式值一定是0

> gaussian elimination

- 從第一個row開始做高斯消去法，只留下主對角線的元素，其他都消為0

```cpp
//判斷是不是那一行的係數全為0
for(j=0;j<n;j++){
    if(abs(augmentMatrix[row][j])>0.0001){
        flag=1;
        break;
    }
}
if(j==n) flag=0;
```

- 因為是宣告成float所以不能用`>0`，而是`>0.0001`的誤差值
- 係數全為0代表此方程組無唯一解-->行列式值為0

> 主對角線元素

```cpp
if(abs(augmentMatrix[row][row])>0.0001){
        flag=1;
        divisor=augmentMatrix[row][row];
    }
```

- 如果主對角線的元素原本就>0.0001那把它變成1就很簡單，直接除上他自己就好了

```cpp
//如果那個算過小趨近於0，那就從下面剩下的行數相對位置中去尋找大於0的元素
else if(abs(augmentMatrix[row][row])<=0.0001){
    swapRow(row);
    swapTimes++;
    if(flag) toOne(row);
}
```
- 紀錄swapTimes的原因是如果行的總交換次數是奇數次那行列式值就是原來的行列式值，如果是奇數次，那行列式值要乘個負號-->美交換一次就乘一個負號
- swapRow是去找下面剩下的row中相對應個元素非0的可以跟這個row進行交換

> 行行互換

```cpp
for(i=row;i<n;i++){
    if(abs(augmentMatrix[i][row])>0.0001){
        for(j=0;j<n+1;j++){
            temp[j]=augmentMatrix[row][j];
            augmentMatrix[row][j]=augmentMatrix[i][j];
            augmentMatrix[i][j]=temp[j];
        }
        flag=1;
        break;
    }
}
```
> 行列式值

```cpp
determinant*=divisor;  //邊做高斯消去邊算行列式值（原本的對角線元素相乘）
for(j=0;j<n+1;j++) augmentMatrix[row][j]/=divisor;
```

- divisor是要把對角線元素變成1時要除的那個數字

> 驗算程式

- 其實很簡單，就是把算出來的向量解代回元方程組Ax=b看等號兩邊是否相等
