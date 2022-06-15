# About this branch

This branch hosts only the file `numerical_result_monad.rds.lrz`. Github doesn't
like commits bigger than 100MB so the file have to be compressed. To decompress
the file you need the `lrzip` [tool](https://github.com/ckolivas/lrzip), which
is usually included in your Linux distribution's 'repository. After installing
`lrzip`, use this command to decompress the file:

```
lrunzip numerical_result_monad.rds.lrz
```

Then use `readRDS` in R to read in the new `.rds` file output by the `lrunzip`
command.
