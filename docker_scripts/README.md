**How to test with a new Trilinos version**
- Pick a Trilinos version to use. You can use a tag or a commit hash.
- Add another image definition to `ci-docker.yml` (see [6231ba9](https://github.com/Pressio/pressio/commit/6231ba9b5d3b11e820951f8cdb83ffcdb0c75577)).
- Enable pushing images from branches to get the new image into registry (see [38a7fb5](https://github.com/Pressio/pressio/commit/38a7fb5e5c03a53aa77b711de53855614582293c)).
  - Normally pushing only happens on `develop`. In this case we need the new image to be pushed from branch to be able to run the tests from `ci-trilinos.yml`.
- You need to let CI run at least once for the image to be pushed to container registry.
- Once the image is published, remember to disable pushing from branches (see [b9540c7](https://github.com/Pressio/pressio/commit/b9540c756bb2c6de3fb24e2c78d3f2c06c524768)).
- You can now enable testing using the newly added image in `ci-trilinos.yml` (see [496df2d](https://github.com/Pressio/pressio/commit/496df2d33f085518e5897baad825c893b588367d)).

See how Trilinos `702aac5` was added in https://github.com/Pressio/pressio/pull/651.

---
**Debugging notes**

You can double-check Trilinos version in `Build and Push Docker images for Github Container Registry` step:
```
#17 46.23 HEAD is now at 702aac58950 Merge Pull Request #11800 from gsjaardema/Trilinos/seacas-snapshot-04-18-23
```
You should also see a successful push at the end of that step:
```
#20 pushing manifest for ghcr.io/pressio/ubuntu-gnu-trilinos-11:702aac5@sha256:0428a42abee7e5c322c0a8527b0653b326bb59486fa7246a90d629742400a330 3.4s done
```
