# Releases

The package name is `IBist`; its GitHub repository is `ibist/IBist.jl`.
The repository is public, but General registration is a separate process.
Until accepted, install by GitHub URL as shown in the README.

## Prepare a version

1. Set `version` in Project.toml to a release such as `0.1.0`, without a
   prerelease suffix. Update CHANGELOG.md for that release.
2. Check dependencies and compatibility bounds, run `Pkg.test()`, and build
   the manual. Keep generated documentation and manifests out of Git.
3. Commit and push. Confirm that package and documentation CI pass.
4. Authorize the JuliaRegistrator GitHub App for this repository if needed:
   https://github.com/apps/juliateam-registrator/installations/new
5. Comment on the release commit:

   ```text
   @JuliaRegistrator register
   ```

Registrator opens a General pull request. Address its checks and any
maintainer feedback. New registrations have a waiting period; package-name
similarity checks may require a maintainer override. Do not infer approval
from a successful registration request.

After General accepts the version, TagBot creates the corresponding tag
and GitHub release. It uses the workflow's GitHub token; no additional secret
is needed for this repository's release-only setup. A maintainer can run
TagBot manually if its normal registry notification was missed.

Official guidance:
- [General registry](https://github.com/JuliaRegistries/General)
- [Registrator](https://github.com/JuliaRegistries/Registrator.jl)
- [TagBot](https://github.com/JuliaRegistries/TagBot)
