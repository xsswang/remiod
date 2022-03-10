# Round 1

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

* This is a new release. No additional notes.

* Reverse dependencies
This is a new release, so there are no reverse dependencies.

# Round 2

## R CMD check results

0 errors | 0 warnings | 0 note

### Reviewer comments
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

Comments are all addressed. 
