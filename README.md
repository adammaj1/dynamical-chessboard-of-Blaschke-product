The dynamical chessboard of finite Blaschke-product

B3(z) = ((z+1/2)/(z/2+1))^3

It is inspired by second version of the paper ["PARABOLIC AND NEAR-PARABOLIC RENORMALIZATION FOR LOCAL DEGREE THREE" by FEI YANG](https://arxiv.org/pdf/1510.00043v2.pdf)


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

