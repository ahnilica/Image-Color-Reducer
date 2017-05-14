# Image-Color-Reducer
A C++ implementation of the winner-takes-all clustering algorithm to approximate the optimal subset of colors found in the image which can be used to reduce its color set. Running with 32 clusters will result in a new image that looks similar to the old image, but only has 32 colors present in it. It uses a matrix library, which has been included. 

## Results

#### original Image:
![Original](/flowers.jpg)

#### 64 Colors:
![64](/examples/output64.jpg)

#### 32 Colors:
![32](/examples/output32.jpg)

#### 16 Colors:
![16](/examples/output16.jpg)

#### 8 Colors:
![8](/examples/output8.jpg)

#### 4 Colors:
![4](/examples/output4.jpg)

#### 2 Colors:
![2](/examples/output2.jpg)

        Usage: NumberOfClusters learningRate filename
        
**note** a good learning rate value is small, such as 0.01


