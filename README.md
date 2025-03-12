# bio-anno-rs

A collection of utilities for reading, writing, and handling records
contained in file formats commonly-used in bioinformatics.

This is very unstable code and should be expected to change
without warning.

## TODO:

* Write files
* Add random shuffling of records (see [random iteration over RecordsDB](https://github.com/jwschroeder3/DNAshape_motif_finder/blob/7f4d46ca9d75b6fbf9006ef651e306e62b37e578/rust_utils/motifer/src/lib.rs#L2671) for start point)
* Add sorting by chromosome name and start position
* Handle more data types.
* Use traits to simplify code for common operations, such as filtering, sorting, and shuffling records in a narrowpeak, bedgraph, gff, etc., file.

## use in a rust binary

To use this crate, I *think* you can add the following to your
`Cargo.toml` file:

```toml
[dependencies]
bio_anno = { git = "https://github.com/jwschroeder3/bio-anno-rs", branch = "main" }
```

And then include `use bio-anno-rs;` in your source code.

