# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given where it's due.

A few things to note before you begin:

- The LASER team is committed to maintaining a welcoming community. We ask that all contributors adhere to our [Code of Conduct](CODE_OF_CONDUCT.md).
- Code contributions should follow the Google Python [Style Guide](https://google.github.io/styleguide/pyguide.html).
- When [reporting a bug](https://github.com/laser-base/laser-core/issues) please include details, such as your operating system and the exact steps to reproduce the issue.
- LASER could always use more documentation, whether as part of the official LASER docs, in docstrings, or even on the web in blog posts, articles, and such.
- The best way to send feedback on LASER is to [file an issue](https://github.com/laser-base/laser-core/issues).

If you have any questions, please reach out to the LASER moderators, [Christopher Lorton](chrisotpher.lorton@gatesfoundation.org) and [Paul Saxman](paul.saxman@gatesfoundation.org).

## Contributing Code

For contriburing code to LASER, you'll need to set up a local development environment:

1. Fork [`laser-generic`](https://github.com/laser-base/laser-generic) (look for the "Fork" button).

2. Create a local clone of your fork:

   ```sh
   git clone git@github.com:YOUR_GITHUB_NAME/laser-generic.git
   ```

3. Create a branch for local development:

   ```sh
   git checkout -b name-of-your-bugfix-or-feature
   ```

4. Make your changes locally.

5. When you're done making changes run all the checks and docs builder with one command:

   ```sh
   tox
   ```

6. Commit your changes and push your branch to GitHub:

   ```sh
   git add .
   git commit -m "Your detailed description of your changes."
   git push origin name-of-your-bugfix-or-feature
   ```

7. Submit a pull request through the GitHub website.

### Pull Request Guidelines

If you need some code review or feedback while you're developing the code just make the pull request.

For merging, you should:

1. Include passing tests (run `tox`).
2. Update documentation when there are new APIs, functionality, etc.
3. Add a note to `CHANGELOG.md` about the changes.
4. Add yourself to `AUTHORS.md`.

### Tips

To run a subset of tests:

```sh
tox -e envname -- pytest -k test_myfeature
```

To run all the test environments in *parallel*:

```sh
tox -p auto
```
