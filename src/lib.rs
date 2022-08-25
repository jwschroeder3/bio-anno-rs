use std::fmt;
use std::io::BufReader;
use std::fs::File;
use std::error::Error;
use std::path;
use serde::Deserialize;
use std::process;
use itertools::Itertools;
use csv;

#[cfg(test)]
mod tests {

    /////////////////////////////////////////////////////////
    // update to work anywhere, setting src depending on where we are?
    /////////////////////////////////////////////////////////
    const TESTDIR: &str = "/home/jeremy/src/bio-anno-rs/test_files";
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::Array;

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
    fn test_get_cpm() {
        let mut bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("cov.bedgraph"),
        ).unwrap();
        let cpm = bgd.get_cpm().unwrap();
        let answer = vec![
            17133.54, 17133.54, 17133.54,
            159880.08, 214033.18, 124028.77,
            307175.93, 99303.01, 44178.41,
        ];
        for (i,res) in cpm.iter().enumerate() {
            assert_abs_diff_eq!(*res, answer[i], epsilon=1e-2);
        }
    }

    #[test]
    fn test_to_cpm() {
        let mut bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("cov.bedgraph"),
        ).unwrap();
        let answer = vec![
            17133.54, 17133.54, 17133.54,
            159880.08, 214033.18, 124028.77,
            307175.93, 99303.01, 44178.41,
        ];
        bgd.to_cpm();
        for (i,res) in answer.iter().enumerate() {
            assert_abs_diff_eq!(*res, bgd.data[i].score, epsilon=1e-2);
        }
    }


    #[test]
    fn test_get_contigs() {
        let bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("test.bedgraph"),
        ).unwrap();
        let answer = vec!["CP064350.1", "CP064351.1", "pBRP02"];
        let contigs = bgd.get_contigs();
        itertools::assert_equal(contigs, answer);
    }

    #[test]
    fn test_res() {
        let bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("test.bedgraph"),
        ).unwrap();
        assert_eq!(5, bgd.get_resolution())
    }

    #[test]
    fn test_ctg_lengths() {
        let bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("test.bedgraph"),
        ).unwrap();
        let contigs = bgd.get_contigs();
        let answer: Vec<usize> = vec![55025, 1070350, 9945];
        let mut results: Vec<usize> = vec![];
        for contig in contigs {
            results.push(bgd.get_contig_length(&contig).unwrap());
        }
        itertools::assert_equal(results, answer);
    }

    #[test]
    fn test_roll_mean() {
        let bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("small.bedgraph"),
        ).unwrap();
        let winsize: usize = 3;
        let answer = vec![
            0.06669717398000229,
            0.06669717398000229,
            0.06669717398000229,
            -0.646127,
            -0.646127,
            -0.646127,
            -0.5847706,
            -0.5847706,
            -0.5847706,
        ];
        let result = bgd.roll_mean(
            winsize,
            true,
        ).unwrap();
        let res_scores = result.fetch_scores().unwrap();
        for (i,res) in res_scores.iter().enumerate() {
            assert_abs_diff_eq!(*res, answer[i], epsilon=1e-5);
        }

        let result = bgd.roll_mean(
            winsize,
            false,
        ).unwrap();
        let res_scores = result.fetch_scores().unwrap();
        let answer = vec![
            0.06669717398000229,
            0.06669717398000229,
            0.06669717398000229,
            -0.6926474,
            -0.646127,
            -0.5996065,
            -0.9260348,
            -0.5847706,
            -0.2435064,
        ];
        for (i,res) in res_scores.iter().enumerate() {
            assert_abs_diff_eq!(*res, answer[i], epsilon=1e-5);
        }
    }

    #[test]
    fn test_padding() {
        let bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("small.bedgraph"),
        ).unwrap();
        let padded = bgd.get_padded_scores(2, true).unwrap();
        let answer = vec![
            -0.386565212080191,
            -0.1719770035449314,
            0.06669717398000229,
            0.06669717398000229,
            0.06669717398000229,
            -0.622378649894243,
            -0.8331850071880074,
            -0.48281724264735265,
            -1.1957696376198523,
            -0.386565212080191,
            -0.1719770035449314,
            0.06669717398000229,
            0.06669717398000229,
        ];
        itertools::assert_equal(padded, answer);

        let padded = bgd.get_padded_scores(2, false).unwrap();
        let answer = vec![
            0.06669717398000229,
            0.06669717398000229,
            0.06669717398000229,
            0.06669717398000229,
            0.06669717398000229,
            -0.622378649894243,
            -0.8331850071880074,
            -0.48281724264735265,
            -1.1957696376198523,
            -0.386565212080191,
            -0.1719770035449314,
            -0.1719770035449314,
            -0.1719770035449314,
        ];
        itertools::assert_equal(padded, answer);
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

    #[test]
    fn test_fetch_scores() {
        let bgd = BEDGraphData::from_file(
            &path::Path::new(TESTDIR).join("test.bedgraph"),
        ).unwrap();
        let scores = bgd.fetch_scores().unwrap();
        let answer = vec![0.06669717398000229; 5];
        assert_eq!(answer, scores[0..5]);
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

    fn set_score(&mut self, new_score: f64) {
        self.score = new_score;
    }
}

/// Implement `Display` for `BEDGraphRecord`.
impl fmt::Display for BEDGraphRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\t{}\t{}\t{}", self.seqname, self.start, self.end, self.score)
    }
}

/// holds a bedgraph file
pub struct BEDGraphData {
    data: Vec<BEDGraphRecord>,
}

impl BEDGraphData {
    /// Read a bedgraph file
    pub fn from_file(fname: &path::PathBuf) -> Result<BEDGraphData, Box<dyn Error>> {

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

    /// Writes the bedgraph file to stdout
    pub fn print(&self) -> Result<(), Box<dyn Error>> {
        for record in &self.data {
            println!("{}", record);
        }
        Ok(())
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

    fn get_cpm(&self) -> Result<Vec<f64>, Box<dyn Error>> {
        let scores = self.fetch_scores()?;
        let sum: f64 = scores.iter().sum();
        let cpm: Vec<f64> = scores.iter()
            .map(|a| a / sum * 1_000_000.0)
            .collect();
        Ok(cpm)
    }

    pub fn to_cpm(&mut self) -> Result<(), Box<dyn Error>> {
        let cpm = self.get_cpm()?;
        for (i,new_score) in cpm.iter().enumerate() {
            self.data[i].set_score(*new_score);
        }
        Ok(())
    }

    /// returns a Result, which if successful, contains the score
    /// column of the bedgraph file as a vector of f64 values
    pub fn fetch_scores(&self) -> Result<Vec<f64>, Box<dyn Error>> {
        let scores: Vec::<f64> = self.iter()
            .map(|record| record.score)
            .collect();
        Ok(scores)
    }

    /// returns a vec of contig names
    fn get_contigs(&self) -> Vec<String> {
        let contigs = self.iter()
            .map(|x| x.seqname)
            .unique()
            .collect();
        contigs
    }

    fn get_max_end(&self) -> usize {
        self.iter().map(|x| x.end).max().unwrap()
    }

    /// returns the greatest end position in self
    fn get_contig_length(&self, seqname: &str) -> Result<usize, Box<dyn Error>> {
        let ctg_bg = self.filter(
            seqname,
            &0,
            &usize::MAX,
        )?;
        Ok(ctg_bg.get_max_end())
    }

    fn get_padded_scores(
            &self,
            pad_size: usize,
            circular: bool,
    ) -> Result<Vec<f64>, Box<dyn Error>> {
        let scores = self.fetch_scores()?;
        let n_scores = scores.len();
        let mut padded: Vec<f64> = Vec::with_capacity(n_scores + pad_size * 2);
        if circular {
            padded.extend_from_slice(&scores[n_scores-pad_size..]);
            padded.extend_from_slice(&scores[..]);
            padded.extend_from_slice(&scores[..pad_size]);
        } else {
            for _ in 0..pad_size {
                padded.push(scores[0]);
            }
            padded.extend_from_slice(&scores[..]);
            for _ in 0..pad_size {
                padded.push(scores[n_scores-1]);
            }
        }
        Ok(padded)
    }

    pub fn get_resolution(&self) -> usize {
        self.data[1].start - self.data[0].start
    }

    /// calculates a rolling mean for each contig in the bedgraph file
    pub fn roll_mean(
            &self,
            window_size: usize,
            circular: bool,
    ) -> Result<BEDGraphData, Box<dyn Error>> {
        if window_size % 2 == 0 {
            eprintln!("Window size should be an odd number for rolling mean, but you entered {}. Exiting now.", window_size);
            process::exit(1);
        }
        let win_size_f = window_size as f64;
        let contigs = self.get_contigs();
        let mut records: Vec::<BEDGraphRecord> = Vec::with_capacity(self.len());
        for contig in contigs {
            let contig_bg = self.filter(
                &contig,
                &0,
                &usize::MAX,
            )?;
            let padded_scores = contig_bg.get_padded_scores(
                (window_size-1)/2,
                circular,
            )?;

            let mut mean_prev_opt: Option<(f64, f64)> = None;
            let mut results: Vec<f64> = Vec::with_capacity(contig_bg.len());

            for window in padded_scores.windows(window_size) {
                let mean = match mean_prev_opt {
                    None => {
                        window.iter().sum::<f64>() / win_size_f
                    },
                    Some((prev_mean, prev)) => {
                        let next = window.last().unwrap();
                        let mean = prev_mean + (*next - prev) / win_size_f;
                        mean
                    }
                };
                let prev = window.first().unwrap();
                results.push(mean);
                mean_prev_opt = Some((mean, *prev))
            }

            for (i,result) in results.iter().enumerate() {
                let record = BEDGraphRecord::new(
                    contig.to_string(),
                    contig_bg.data[i].start,
                    contig_bg.data[i].end,
                    *result,
                );
                records.push( record );
            }
        }
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

