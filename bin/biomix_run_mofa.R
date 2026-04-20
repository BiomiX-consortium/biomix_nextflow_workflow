cli_args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag) {
  position <- match(flag, cli_args)
  if (is.na(position) || position == length(cli_args)) {
    stop(sprintf("Missing required argument: %s", flag), call. = FALSE)
  }
  cli_args[[position + 1]]
}

workspace <- normalizePath(parse_arg("--workspace"), winslash = "/", mustWork = TRUE)
group1 <- parse_arg("--group1")
group2 <- parse_arg("--group2")

setwd(workspace)
combined_json <- jsonlite::fromJSON(
  txt = paste(readLines(file.path(workspace, "COMBINED_COMMANDS.json"), warn = FALSE), collapse = "\n")
)

assign("args", as.list(c(group1, group2, workspace)), envir = .GlobalEnv)
assign("directory", workspace, envir = .GlobalEnv)
assign("combined_json", combined_json, envir = .GlobalEnv)
assign("COMMAND", combined_json[["COMMANDS"]], envir = .GlobalEnv)
assign("COMMAND_MOFA", combined_json[["COMMANDS_MOFA"]], envir = .GlobalEnv)
assign("COMMAND_ADVANCED", combined_json[["COMMANDS_ADVANCED"]], envir = .GlobalEnv)
assign("DIR_METADATA", combined_json[["DIRECTORY_INFO"]][["METADATA_DIR"]], envir = .GlobalEnv)
assign("DIR_METADATA_output", combined_json[["DIRECTORY_INFO"]][["OUTPUT_DIR"]], envir = .GlobalEnv)

source(file.path(workspace, "Integration", "MOFA_MULTI2.R"), local = globalenv())
