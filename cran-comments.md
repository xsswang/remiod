# remiod (version 1.0.1)

## Round 1

### R CMD check results
0 errors | 0 warnings | 2 notes

NOTEs:

Found the following (possibly) invalid URLs:
  URL: https://cran.r-project.org/web/packages/remiod/index.html
    From: README.md
    Status: 200
    Message: OK
    CRAN URL not in canonical form
  URL: https://github.com/xsswang/remiod/main/LICENSE
    From: README.md
    Status: 404
    Message: Not Found
  The canonical URL of the CRAN page for a package is 
    https://CRAN.R-project.org/package=pkgname


Found if() conditions comparing class() to string:
File 'remiod/R/get_modeltypes_custom.R': if (!is.null(auxvars) & class(auxvars) != "formula") ...
File 'remiod/R/model_imp_custom.R': if (n.iter > 0 & class(mcmc) != "mcmc.list") ...

REPLY:

URLs are revised.

for comparing class() to string, it is replaced with inherits().


# remiod (version 1.0.0)

## Round 1

### R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

* This is a new release. No additional notes.

* Reverse dependencies
This is a new release, so there are no reverse dependencies.

## Round 2

### R CMD check results

0 errors | 0 warnings | 0 note

#### Reviewer comments
2022-03-09 Julia Haider

```
Thanks, in the description text you have

framework.Methods

Please add a space:

framework. Methods

If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form authors (year) <doi:...> authors (year) <arXiv:...> authors (year, ISBN:...) or if those are not available: <https:...> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

Please do not modify the global environment in your functions. This is not allowed by the CRAN policies.
.Random.seed should not be altered by the user.
```


2022-03-10 Gregor Seyer

```
Thanks, in the documentation. Please write about the structure of the output
(class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar) Missing Rd-tags:
      get_Mlist.Rd: \value
      list.models.Rd: \value
      mcmcplot.Rd: \value
      summary.Rd: \value

Please add a few more small executable examples in your Rd-files to illustrate the use of the exported function but also enable automatic testing.
```

Comments are all addressed. 

