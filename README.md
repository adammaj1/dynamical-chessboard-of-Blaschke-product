The dynamical chessboard of finite Blaschke-product

B3(z) = ((z+1/2)/(z/2+1))^3

It is inspired by second version of the paper ["PARABOLIC AND NEAR-PARABOLIC RENORMALIZATION FOR LOCAL DEGREE THREE" by FEI YANG](https://arxiv.org/pdf/1510.00043v2.pdf),  

Here is original image extracted from the paper, see Figure 14 from page 52.: " the first column indicate the dynamical chessboards of f and B3" 

![](./png/35.png "original image") 


Here is my image, it looks like 1/3 of original image:

![](./png/IBD_LCM.png "my image") 


# Analysis
* B3 is a Blaschke product whose Julia set is the unit circle
* The point z = 1 is a 1-parabolic fixed point with two attracting petals
* the unit disk D and C∖D are two immediate basins of 1
* critical point:  only one point z=-1/2
* poles: z = -2 (order 3 pole)
* zeros: z = -1/2 (order 3 zero)


# Discussion
* [math.stackexchange question: julia-set-of-finite-blaschke-product](https://math.stackexchange.com/questions/4448405/julia-set-of-finite-blaschke-product)
* [fractalforums.org: critical-points-of](https://fractalforums.org/fractal-mathematics-and-new-theories/28/critical-points-of/4753)


# Algorithm


# Source code
Program: 
* [j.c](./src/j.c ) - c console program
* [j.sh](./src/j.sh ) - bash script
* [Makefile](./src/Makefile) - Makefile

Run: simply make from console


Output: 
* [png images](./png/) 
* [j.txt](./src/j.txt ) - text output



# other images

![](./png/FatouBasins_LCM.png)  

![](./png/FatouBasins_LSCM.png)  

![](./png/FatouBasins_LSCM_trap.png)  
![](./png/IBD.png)  

![](./png/IBD_LSCM.png)  

![](./png/IBD_LSM_LCM.png)  

![](./png/IBD_LSM_LSCM.png)  

![](./png/LCM.png)  

![](./png/LS2CM.png)  

![](./png/LS2CM_cr.png)
   
![](./png/LS2M.png)  

![](./png/LSCM.png)  

![](./png/LSM.png)  

![](./png/MBD.png)  

![](./png/MBD_LCM.png  )  

![](./png/MBD_LSCM.png)   

![](./png/MBD_LSM_LCM.png) 

![](./png/MBD_LSM_LSCM.png) 

![](./png/ParabolicCheckerboard2_LCM.png) 

![](./png/ParabolicCheckerboard2_LSCM2.png) 

![](./png/ParabolicCheckerboard2_LSM.png) 

![](./png/ParabolicCheckerboard_LSCM.png)  

![](./png/ParabolicCheckerboard_LSM.png)





# Git

create a new repository on the command line

```git
echo "# The dynamical chessboard of finite Blaschke-product" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:adammaj1/dynamical-chessboard-of-Blaschke-product.git
git push -u origin main
```


## Local repo
```
~/Dokumenty/julia_blaszke/b3/dynamical-chessboard-of-Blaschke-product/
```




## Subdirectory

```git
mkdir png
git add *.png
git mv  *.png ./png
git commit -m "move"
git push -u origin main
```
then link the images:

```txt
![](./png/n.png "description") 

```

```git
gitm mv -f 
```

[Remove a file/directory from a Git repository without deleting it from the local filesystem](https://stackoverflow.com/questions/1143796/remove-a-file-from-a-git-repository-without-deleting-it-from-the-local-filesyste)

```git
git rm --cached myfile.log
```

Single directory and all it's content

```git
git rm --cached -r ./png
git commit -m "png"
git push -u origin main

```


Rename directory

```git
git add ./lcm/
git mv ./lcm ./png 
git commit -m "png"
git push -u origin main

```

## Github
* [GitHub Flavored Markdown Spec](https://github.github.com/gfm/)
* [md cheat sheet](http://mdcheatsheet.com/)
* [CommonMark Spec](https://spec.commonmark.org)
* [Markdown parser ](https://markdown-it.github.io/)

