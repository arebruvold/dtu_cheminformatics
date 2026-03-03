# ms scripts
library(tidyverse)


pcprop_from_name_vec <-
  function(
    in_name,
    property = c(
      "XLogP",
      "ConnectivitySMILES",
      "MolecularFormula",
      "MolecularWeight",
      "ExactMass",
      "TPSA",
      "HeavyAtomCount",
      "Charge",
      "InChI",
      "InChIKey",
      "IUPACName"
    ),
    cache_file = "pubchem_all_prop_cache.rds"
  ) {
    # --- 1. Define ALL Properties and Desired Property List ---

    # The complete list of properties to ALWAYS request from PubChem
    # NOTE: Explicitly including "CID" here ensures it is handled in filtering,
    # and "XLogP" which was missing from the default argument but present in ALL_PROPERTIES.
    ALL_PROPERTIES <- c(
      "CID", # <-- ADDED: Explicitly handle the implicitly returned CID
      "XLogP",
      "ConnectivitySMILES",
      "MolecularFormula",
      "MolecularWeight",
      "ExactMass",
      "TPSA",
      "HeavyAtomCount",
      "Charge",
      "InChI",
      "InChIKey",
      "IUPACName"
    )

    # Standardize the user's requested property list for output selection
    if (length(property) == 1 && property == "all") {
      requested_properties <- ALL_PROPERTIES
    } else {
      # Ensure requested properties are valid and available in ALL_PROPERTIES
      invalid_props <- setdiff(property, ALL_PROPERTIES)
      if (length(invalid_props) > 0) {
        warning(
          sprintf(
            "Ignoring invalid properties: %s",
            paste(invalid_props, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      requested_properties <- intersect(property, ALL_PROPERTIES)
      if (length(requested_properties) == 0) {
        stop(
          "The 'property' argument must contain at least one valid property.",
          call. = FALSE
        )
      }
    }

    # The property string used for the API call (always ALL properties *excluding* CID,
    # since CID is implicit, but let's just use the existing logic for safety)
    # The actual API string must NOT contain "CID" if we want to retrieve the properties correctly
    API_PROPS_FOR_CURL <- setdiff(ALL_PROPERTIES, "CID")
    API_PROPERTY_STRING <- stringr::str_c(API_PROPS_FOR_CURL, collapse = ",")

    # --- 2. Load or Initialize Cache (No change) ---
    if (file.exists(cache_file)) {
      local_cache <- tryCatch(
        readRDS(cache_file),
        error = function(e) {
          warning(
            sprintf(
              "Failed to read cache file: %s. Starting empty.",
              conditionMessage(e)
            ),
            call. = FALSE
          )
          return(data.frame(
            name = character(0),
            result = list(),
            stringsAsFactors = FALSE
          ))
        }
      )
    } else {
      message(sprintf("Cache file '%s' not found. Starting empty.", cache_file))
      local_cache <- data.frame(
        name = character(0),
        result = list(),
        stringsAsFactors = FALSE
      )
    }

    # Check for valid cache structure
    if (
      !is.data.frame(local_cache) ||
        !all(c("name", "result") %in% names(local_cache)) ||
        !is.list(local_cache$result)
    ) {
      warning("Cache structure invalid. Resetting cache.", call. = FALSE)
      local_cache <- data.frame(
        name = character(0),
        result = list(),
        stringsAsFactors = FALSE
      )
    }

    # --- 3. Internal Worker Function (Cache-aware, always fetches ALL) ---
    .single_pcprop_from_name_cached <- function(single_name, current_cache) {
      # 3.1 Check Cache (No change)
      cached_row_index <- which(current_cache$name == single_name)

      if (length(cached_row_index) > 0) {
        # Cache hit: Return the ALL properties result
        return(list(
          result_all = current_cache$result[[cached_row_index[1]]],
          new_entry = data.frame() # No new entry needed
        ))
      }

      # 3.2 EXTERNAL API LOOKUP
      escaped_name <- stringr::str_replace_all(single_name, " ", "%20")
      out_type <- "/CSV"

      # Construct and execute the curl command
      out <- system(
        paste0(
          'curl -g "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/',
          escaped_name,
          "/property/",
          API_PROPERTY_STRING,
          out_type,
          '"'
        ),
        intern = TRUE,
        ignore.stderr = TRUE
      )

      # Process output first (before handling multiple matches)
      out_df <- tryCatch(
        {
          read.csv(text = out, stringsAsFactors = FALSE)
        },
        error = function(e) data.frame()
      )

      # Handle multiple matches warning ONLY IF we have valid data to cache
      if (length(out) > 2) {
        # Check if we will actually cache this result
        will_cache <- FALSE

        if (nrow(out_df) > 0) {
          property_cols <- intersect(ALL_PROPERTIES, names(out_df))
          if (length(property_cols) > 0) {
            value_cols <- setdiff(property_cols, "CID")
            if (length(value_cols) > 0) {
              all_na <- all(is.na(out_df[, value_cols, drop = FALSE]))
              will_cache <- !all_na
            }
          }
        }

        # Only warn about multiple matches if we have valid data to cache
        if (will_cache) {
          warning(
            paste(
              "Warning! Multiple matches for cname:",
              single_name,
              ". Using first."
            ),
            call. = FALSE
          )
        }

        # Always truncate to first result regardless
        out <- out[1:2]

        # Re-process the truncated output
        out_df <- tryCatch(
          {
            read.csv(text = out, stringsAsFactors = FALSE)
          },
          error = function(e) data.frame()
        )
      }

      # Respect API limit
      Sys.sleep(0.25)

      # Process output
      out_df <- tryCatch(
        {
          read.csv(text = out, stringsAsFactors = FALSE)
        },
        error = function(e) data.frame()
      )

      # ----------------------------------------------------
      # 3.3 Cache Update Logic - ONLY CACHE IF DATA IS VALID
      # ----------------------------------------------------
      new_entry <- data.frame()

      # Check if the result contains valid data
      if (nrow(out_df) > 0) {
        # Check for valid columns that came back, including the implicit CID
        # property_cols now includes CID, which is necessary for the next step.
        property_cols <- intersect(ALL_PROPERTIES, names(out_df))

        if (length(property_cols) > 0) {
          # Check if all actual *property* values (excluding CID) are NA
          # Filter property_cols to exclude CID for the all_na check
          value_cols <- setdiff(property_cols, "CID")

          if (length(value_cols) > 0) {
            all_na <- all(is.na(out_df[, value_cols, drop = FALSE]))
          } else {
            # Only CID was returned, which isn't a property we check for NA value
            all_na <- FALSE
          }

          if (all_na) {
            warning(
              sprintf(
                "❌ No valid properties found for compound '%s' - all values are NA. Not caching.",
                single_name
              ),
              call. = FALSE
            )
          } else {
            # Cache the result (which implicitly includes CID)
            new_entry <- data.frame(
              name = single_name,
              stringsAsFactors = FALSE
            )
            new_entry$result <- list(out_df)
            message(sprintf("Caching valid result for '%s'", single_name))
          }
        } else {
          warning(
            sprintf(
              "❌ No recognizable property columns found for compound '%s'. API response may be malformed.",
              single_name
            ),
            call. = FALSE
          )
        }
      } else {
        warning(
          sprintf("❌ No data returned for compound '%s'.", single_name),
          call. = FALSE
        )
      }

      return(list(result_all = out_df, new_entry = new_entry))
    }

    # --- 4. Main Caching and Vectorization Logic (No change in logic) ---
    # ... (Omitted for brevity, as the logic remains the same, but relies on changes in Step 3) ...

    if (length(in_name) == 0) {
      return(NULL)
    }
    unique_names <- unique(in_name)
    n_total <- length(unique_names)
    new_cache_entries <- data.frame(
      name = character(0),
      result = list(),
      stringsAsFactors = FALSE
    )
    all_prop_results_list <- vector("list", n_total)

    message(sprintf(
      "Starting lookup for %d unique compounds. Cache file: %s",
      n_total,
      cache_file
    ))

    for (i in seq_along(unique_names)) {
      single_name <- unique_names[i]
      res <- .single_pcprop_from_name_cached(single_name, local_cache)
      all_prop_results_list[[i]] <- res$result_all

      if (nrow(res$new_entry) > 0) {
        new_cache_entries <- rbind(new_cache_entries, res$new_entry)
        message(sprintf("Added pc properties %d/%d", i, n_total))
      } else if (nrow(res$result_all) > 0) {
        message(sprintf("Used cached properties %d/%d", i, n_total))
      } else {
        message(sprintf(
          "Failed lookup %d/%d for '%s'",
          i,
          n_total,
          single_name
        ))
      }
    }

    # 5. Final Cache Update and Save (No change in logic)
    if (nrow(new_cache_entries) > 0) {
      current_names <- local_cache$name
      new_entries_unique <- new_cache_entries[
        !(new_cache_entries$name %in% current_names),
      ]

      if (nrow(new_entries_unique) > 0) {
        updated_cache <- rbind(local_cache, new_entries_unique)
        saveRDS(updated_cache, file = cache_file)
        message(sprintf(
          "✅ Saved %d new valid query results to cache file '%s'.",
          nrow(new_entries_unique),
          cache_file
        ))
      }
    }

    # --- 6. Format Final Output (Map results back to original input vector) ---

    # Create a lookup table from unique names to their full property results
    name_to_all_props <- all_prop_results_list
    names(name_to_all_props) <- unique_names
    final_output_list <- name_to_all_props[in_name]

    # CID is implicitly handled in the output by adding it to requested_properties
    final_requested_properties <- unique(c("CID", requested_properties))

    # Determine if the output should be a single vector or a data frame
    if (length(requested_properties) == 1) {
      # Single property requested -> return a vector (CID is not forced in single output)
      output_vec <- purrr::map_chr(
        final_output_list,
        ~ {
          if (
            is.data.frame(.) &&
              nrow(.) > 0 &&
              requested_properties %in% names(.)
          ) {
            return(as.character(.[[requested_properties]][1]))
          } else {
            return(NA_character_)
          }
        }
      )
      return(output_vec)
    } else {
      # Multiple properties requested -> return a combined data frame (CID is forced here)
      output_df <- purrr::map_dfr(
        final_output_list,
        ~ {
          if (is.data.frame(.) && nrow(.) > 0) {
            # Select only the columns present in the result, including the forced CID
            cols_to_select <- intersect(final_requested_properties, names(.))
            df_out <- dplyr::select(., tidyselect::all_of(cols_to_select))

            # Add back missing required columns with NA for consistency
            missing_cols <- setdiff(final_requested_properties, names(df_out))
            for (mc in missing_cols) {
              df_out[[mc]] <- NA
            }
            return(df_out)
          } else {
            # Return an empty row with the required columns set to NA
            df_empty <- setNames(
              data.frame(matrix(
                ncol = length(final_requested_properties),
                nrow = 1
              )),
              final_requested_properties
            )
            df_empty[] <- NA
            return(df_empty)
          }
        }
      )
      return(output_df)
    }
  }


pcprop_from_inchikey_vec <-
  function(
    in_inchikey,
    property = c(
      "XLogP",
      "CanonicalSMILES",
      "MolecularFormula",
      "MolecularWeight",
      "ExactMass",
      "TPSA",
      "HeavyAtomCount",
      "Charge",
      "InChI",
      "InChIKey",
      "IUPACName"
    ),
    cache_file = "pubchem_inchikey_prop_cache.rds"
  ) {
    # --- 1. Define ALL Properties and Desired Property List ---

    # The complete list of properties to ALWAYS request from PubChem
    ALL_PROPERTIES <- c(
      "CID", # Implicitly returned
      "XLogP",
      "CanonicalSMILES",
      "MolecularFormula",
      "MolecularWeight",
      "ExactMass",
      "TPSA",
      "HeavyAtomCount",
      "Charge",
      "InChI",
      "InChIKey",
      "IUPACName"
    )

    # Standardize the user's requested property list
    if (length(property) == 1 && property == "all") {
      requested_properties <- ALL_PROPERTIES
    } else {
      invalid_props <- setdiff(property, ALL_PROPERTIES)
      if (length(invalid_props) > 0) {
        warning(
          sprintf(
            "Ignoring invalid properties: %s",
            paste(invalid_props, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      requested_properties <- intersect(property, ALL_PROPERTIES)
      if (length(requested_properties) == 0) {
        stop(
          "The 'property' argument must contain at least one valid property.",
          call. = FALSE
        )
      }
    }

    # API string excludes CID (it's implicit)
    API_PROPS_FOR_CURL <- setdiff(ALL_PROPERTIES, "CID")
    API_PROPERTY_STRING <- stringr::str_c(API_PROPS_FOR_CURL, collapse = ",")

    # --- 2. Load or Initialize Cache ---
    if (file.exists(cache_file)) {
      local_cache <- tryCatch(
        readRDS(cache_file),
        error = function(e) {
          warning(
            sprintf(
              "Failed to read cache file: %s. Starting empty.",
              conditionMessage(e)
            ),
            call. = FALSE
          )
          return(data.frame(
            name = character(0),
            result = list(),
            stringsAsFactors = FALSE
          ))
        }
      )
    } else {
      message(sprintf("Cache file '%s' not found. Starting empty.", cache_file))
      local_cache <- data.frame(
        name = character(0),
        result = list(),
        stringsAsFactors = FALSE
      )
    }

    if (
      !is.data.frame(local_cache) ||
        !all(c("name", "result") %in% names(local_cache))
    ) {
      local_cache <- data.frame(
        name = character(0),
        result = list(),
        stringsAsFactors = FALSE
      )
    }

    # --- 3. Internal Worker Function (InChIKey-specific) ---
    .single_pcprop_from_inchikey_cached <- function(
      single_inchikey,
      current_cache
    ) {
      # 3.1 Check Cache
      cached_row_index <- which(current_cache$name == single_inchikey)

      if (length(cached_row_index) > 0) {
        return(list(
          result_all = current_cache$result[[cached_row_index[1]]],
          new_entry = data.frame()
        ))
      }

      # 3.2 EXTERNAL API LOOKUP
      # Sanitize input: InChIKeys should not contain spaces, but trim just in case.
      # They are safe for URLs otherwise (only A-Z and hyphens).
      escaped_key <- trimws(single_inchikey)
      out_type <- "/CSV"

      # Construct command using the /compound/inchikey/ endpoint
      out <- system(
        paste0(
          'curl -g "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/',
          escaped_key,
          "/property/",
          API_PROPERTY_STRING,
          out_type,
          '"'
        ),
        intern = TRUE,
        ignore.stderr = TRUE
      )

      # Process output
      out_df <- tryCatch(
        {
          read.csv(text = out, stringsAsFactors = FALSE)
        },
        error = function(e) data.frame()
      )

      # Handle multiple matches (Rare for InChIKey as it's a hash, but good practice)
      if (length(out) > 2) {
        will_cache <- FALSE
        if (nrow(out_df) > 0) {
          property_cols <- intersect(ALL_PROPERTIES, names(out_df))
          if (length(setdiff(property_cols, "CID")) > 0) {
            will_cache <- !all(is.na(out_df[,
              setdiff(property_cols, "CID"),
              drop = FALSE
            ]))
          }
        }

        if (will_cache) {
          warning(
            paste(
              "Warning! Multiple matches for InChIKey. Using first result."
            ),
            call. = FALSE
          )
        }

        out <- out[1:2] # Truncate to header + first row
        out_df <- tryCatch(
          {
            read.csv(text = out, stringsAsFactors = FALSE)
          },
          error = function(e) data.frame()
        )
      }

      # Respect API limit
      Sys.sleep(0.25)

      # 3.3 Cache Update Logic
      new_entry <- data.frame()

      if (nrow(out_df) > 0) {
        property_cols <- intersect(ALL_PROPERTIES, names(out_df))
        if (length(property_cols) > 0) {
          value_cols <- setdiff(property_cols, "CID")

          # Check for validity (not all NA)
          all_na <- if (length(value_cols) > 0) {
            all(is.na(out_df[, value_cols, drop = FALSE]))
          } else {
            FALSE
          }

          if (all_na) {
            # No valid data found
          } else {
            # Cache the result
            new_entry <- data.frame(
              name = single_inchikey,
              stringsAsFactors = FALSE
            )
            new_entry$result <- list(out_df)
            message(sprintf("Caching valid result"))
          }
        }
      }

      return(list(result_all = out_df, new_entry = new_entry))
    }

    # --- 4. Main Caching and Vectorization Logic ---

    if (length(in_inchikey) == 0) {
      return(NULL)
    }
    unique_keys <- unique(in_inchikey)
    n_total <- length(unique_keys)
    new_cache_entries <- data.frame(
      name = character(0),
      result = list(),
      stringsAsFactors = FALSE
    )
    all_prop_results_list <- vector("list", n_total)

    message(sprintf(
      "Starting lookup for %d unique InChIKeys. Cache file: %s",
      n_total,
      cache_file
    ))

    for (i in seq_along(unique_keys)) {
      single_key <- unique_keys[i]
      res <- .single_pcprop_from_inchikey_cached(single_key, local_cache)
      all_prop_results_list[[i]] <- res$result_all

      if (nrow(res$new_entry) > 0) {
        new_cache_entries <- rbind(new_cache_entries, res$new_entry)
        message(sprintf("Added properties %d/%d", i, n_total))
      } else if (nrow(res$result_all) > 0) {
        message(sprintf("Used cached properties %d/%d", i, n_total))
      } else {
        message(sprintf("Failed lookup %d/%d", i, n_total))
      }
    }

    # 5. Final Cache Update and Save
    if (nrow(new_cache_entries) > 0) {
      current_names <- local_cache$name
      new_entries_unique <- new_cache_entries[
        !(new_cache_entries$name %in% current_names),
      ]

      if (nrow(new_entries_unique) > 0) {
        updated_cache <- rbind(local_cache, new_entries_unique)
        saveRDS(updated_cache, file = cache_file)
        message(sprintf(
          "✅ Saved %d new valid results to cache.",
          nrow(new_entries_unique)
        ))
      }
    }

    # --- 6. Format Final Output ---

    name_to_all_props <- all_prop_results_list
    names(name_to_all_props) <- unique_keys
    final_output_list <- name_to_all_props[in_inchikey]

    final_requested_properties <- unique(c("CID", requested_properties))

    if (length(requested_properties) == 1) {
      output_vec <- purrr::map_chr(
        final_output_list,
        ~ {
          if (
            is.data.frame(.) &&
              nrow(.) > 0 &&
              requested_properties %in% names(.)
          ) {
            return(as.character(.[[requested_properties]][1]))
          } else {
            return(NA_character_)
          }
        }
      )
      return(output_vec)
    } else {
      output_df <- purrr::map_dfr(
        final_output_list,
        ~ {
          if (is.data.frame(.) && nrow(.) > 0) {
            cols_to_select <- intersect(final_requested_properties, names(.))
            df_out <- dplyr::select(., tidyselect::all_of(cols_to_select))

            missing_cols <- setdiff(final_requested_properties, names(df_out))
            for (mc in missing_cols) {
              df_out[[mc]] <- NA
            }
            return(df_out)
          } else {
            df_empty <- setNames(
              data.frame(matrix(
                ncol = length(final_requested_properties),
                nrow = 1
              )),
              final_requested_properties
            )
            df_empty[] <- NA
            return(df_empty)
          }
        }
      )
      return(output_df)
    }
  }
