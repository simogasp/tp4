1. Copy the file `post-xxx-sample.txt `(which is in the same folder of
your TEX distribution as this pdf) into the git hooks directory in your
working copy. In our example case, you should end up with a file called
`~/compsci/.git/hooks/post-checkout`

2. If you’re using a unix-like system, don’t forget to make the file executable. 
   Just how you do this is outside the scope of this manual, but one possible way is with commands such as this:
    ```
    chmod g+x post-checkout.
    ```

3. Test your setup with `git checkout master` (or another suitable branch
name). This should generate copies of `gitHeadInfo.gin` in the directories you intended.

4. Now make two more copies of this file in the same directory (hooks),
calling them `post-commit` and `post-merge`, and you’re done. 
As before, users of unix-like systems should ensure these files are marked as executable.