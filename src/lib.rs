use std::io::BufReader;
use std::fs::File;
use std::error::Error;
use std::path;
use serde::Deserialize;
use std::process;
use csv;

#[cfg(test)]
mod tests {

    const TESTDIR: &str = "/home/jeremy/src/bio-anno-rs/test_files";
    use super::*;

    #[test]
    fn test_bg_filter() {
        let bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("test.bedgraph"),
        ).unwrap();
        let filter_bgd = bgd.filter(
            "CP064350.1",
            &50000,
            &55000,
        ).unwrap();
        let answer1 = BEDGraphRecord::new(
            "CP064350.1".to_string(),
            50000,
            50005,
            -0.07215033236083573,
        );
        println!("First record in file: {:?}", bgd[0]);
        println!("First record in filtered data: {:?}", filter_bgd[0]);
        assert_eq!(filter_bgd[0], answer1);

        let filter_bgd = bgd.filter(
            "pBRP02",
            &9802,
            &9883,
        ).unwrap();
        let answer2 = BEDGraphRecord::new(
            "pBRP02".to_string(),
            9805,
            9810,
            -0.21729724583100157,
        );
        assert_eq!(filter_bgd[0], answer2);

        let answer3 = BEDGraphRecord::new(
            "pBRP02".to_string(),
            9875,
            9880,
            0.5703630173280463,
        );
        assert_eq!(filter_bgd[filter_bgd.len()-1], answer3);
    }

    #[test]
    fn test_read_bedgraph() {
        let bgd = BEDGraphData::from_file(
                &path::Path::new(TESTDIR).join("test.bedgraph"),
            ).unwrap();
        let answer1 = BEDGraphRecord::new(
            "CP064350.1".to_string(),
            0,
            5,
            0.06669717398000229,
        );
        assert_eq!(bgd[0], answer1);

        let answer2 = BEDGraphRecord::new(
            "CP064350.1".to_string(),
            25000,
            25005,
            0.25470646448672307,
        );
        assert_eq!(bgd[5000], answer2);
    }
}

/// struct to define a single line of a bedgraph file
#[derive(Debug, Deserialize, PartialEq)]
pub struct BEDGraphRecord {
    seqname: String,
    start: usize,
    end: usize,
    score: f64,
}

impl BEDGraphRecord {
    fn new(
        seqname: String,
        start: usize,
        end: usize,
        score: f64,
    ) -> BEDGraphRecord {
        BEDGraphRecord {
            seqname,
            start,
            end,
            score,
        }
    }
}

/// holds a bedgraph file
pub struct BEDGraphData {
    data: Vec<BEDGraphRecord>,
}

impl BEDGraphData {
    /// Read a bedgraph file
    fn from_file(fname: &path::PathBuf) -> Result<BEDGraphData, Box<dyn Error>> {

        let file = File::open(fname).unwrap_or_else(|err| {
            eprintln!("Problem reading bedgraph file {:?}: {}", fname, err);
            process::exit(1);
        });
        // open buffered reader to bedgraph file
        let buf_reader = BufReader::new(file);

        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(buf_reader);
        let mut records: Vec<BEDGraphRecord> = Vec::new();

        for result in rdr.deserialize() {
            let record = result.unwrap_or_else(|err| {
                eprintln!("Problem with your bedgraph records. Is {:?} a properly-formed bedgraph file?: {}", fname, err);
                process::exit(1);
            });
            records.push(record);
        }
        Ok(BEDGraphData{ data: records })
    }

    /// Returns number of records in self
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Iterate over each record in the database as a [StrandedSequence] value pair
    fn iter(&self) -> BEDGraphDataIter{
        BEDGraphDataIter { loc: 0, bgd: &self, size: self.len() }
    }
    
    /// filters records in self to those within the given range
    /// returns a new BEDGraphData instance
    pub fn filter(
            &self,
            seqname: &str,
            start: &usize,
            end: &usize,
    ) -> Result<BEDGraphData, Box<dyn Error>> {

        let records: Vec::<BEDGraphRecord> = self.iter()
            .filter(|x| {
                x.seqname == seqname
                && x.start >= *start
                && x.end <= *end
            }).collect();
        
        Ok(BEDGraphData{data: records})
    }
}

/// enables slicing of BEDGraphData struct
impl<Idx> std::ops::Index<Idx> for BEDGraphData
where
    Idx: std::slice::SliceIndex<[BEDGraphRecord]>,
{
    type Output = Idx::Output;

    fn index(&self, index: Idx) -> &Self::Output {
        &self.data[index]
    }
}

/// Allows for iteration over a BEDGraphData struct
///
/// # Fields
///
/// * `loc` - Current location in the database
/// * `bgd` - A reference to the [BEDGraphData]
/// * `size` - The number of BEDGraphRecords in bgd
pub struct BEDGraphDataIter<'a> {
    loc: usize,
    bgd: &'a BEDGraphData,
    size: usize,
}

/// Enables iteration over the BEDGraphData. Returns a [BEDGraphRecord] as 
/// each item.
impl<'a> Iterator for BEDGraphDataIter<'a> {
    type Item = BEDGraphRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if self.loc == self.size{
            None
        } else {
            let out_rec = &self.bgd[self.loc];
            self.loc += 1;
            Some(BEDGraphRecord::new(
                out_rec.seqname.to_string(),
                out_rec.start,
                out_rec.end,
                out_rec.score,
            ))
        }
    }
}

