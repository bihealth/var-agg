/// Helper code for BCF/VCF files.
use slog::Logger;

use rust_htslib;
use std::ffi;

use super::errors::*;

/// Build index file for the VCF/BCF file at `path`.
pub fn build_index(logger: &Logger, path: &String) -> Result<()> {
    if path.ends_with(".vcf.gz") {
        info!(logger, "Writing .tbi index file...");
        // TODO: the unsafe stuff should go into rust_htslib
        use rust_htslib;
        use std::ffi;
        let res = unsafe {
            rust_htslib::htslib::tbx_index_build(
                ffi::CString::new(path.as_bytes()).unwrap().as_ptr(),
                0,
                &rust_htslib::htslib::tbx_conf_vcf,
            )
        };
        if res != 0 {
            bail!("Could not create .tbi index");
        }
    } else if path.ends_with(".bcf") {
        info!(logger, "Writing .csi index file...");
        let res = unsafe {
            rust_htslib::htslib::bcf_index_build(
                ffi::CString::new(path.as_bytes()).unwrap().as_ptr(),
                14,
            )
        };
        if res != 0 {
            bail!("Could not create .csi index");
        }
    } else {
        info!(
            logger,
            "Not building index, output file does not having ending .bcf or .vcf.gz"
        );
    }

    Ok(())
}
