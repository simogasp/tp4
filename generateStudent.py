#!/usr/bin/env python3

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s"
)
logger = logging.getLogger()


def main(working_dir: str, archive_name: str, skip_cleaning: bool = False, skip_packaging: bool = False):

    if not os.path.exists(working_dir):
        logger.info("Creating working directory %s", working_dir)
        os.makedirs(working_dir)

    dest_name = os.path.join(working_dir, archive_name)

    logger.info("Creating directory %s", dest_name)
    os.mkdir(dest_name)

    files_to_copy = ["BUILD.md", "CMakeLists.txt", "README.md", "LICENSE"]
    logger.info("Copying files to %s", dest_name)
    for file in files_to_copy:
        shutil.copy(file, dest_name)

    logger.info("Copying directories data and freeglut to %s", dest_name)
    shutil.copytree("src", os.path.join(dest_name, "src"))
    shutil.copytree("data", os.path.join(dest_name, "data"))
    shutil.copytree("freeglut", os.path.join(dest_name, "freeglut"))

    studentify_dir = os.path.join(tempfile.gettempdir(), "tpt")
    logger.info("Cloning studentify.py to %s", studentify_dir)
    subprocess.check_call(["git", "clone", "https://github.com/simogasp/studentipy.git", studentify_dir])

    files_to_studentify = ["src/ObjModel.cpp", "src/rendering.cpp", "src/geometry.cpp", "src/objReader.cpp"]
    for file in files_to_studentify:
        file_cpp = os.path.join(dest_name, file)
        logger.info("Applying studentify to %s", file_cpp)
        subprocess.check_call(
            ["python3", os.path.join(studentify_dir, "studentify.py"), file_cpp, "-o", file_cpp, "--force"])

    logger.info("Removing studentify.py")
    shutil.rmtree(studentify_dir)

    if not skip_packaging:
        logger.info("Generating archive %s.zip in %s", archive_name, working_dir)
        with zipfile.ZipFile(os.path.join(working_dir, f"{archive_name}.zip"), "w") as zip_file:
            for root, dirs, files in os.walk(dest_name):
                for file in files:
                    zip_file.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), dest_name))

    if skip_cleaning:
        return
    logger.info("Cleaning up %s", dest_name)
    shutil.rmtree(dest_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script used to generate the package to be given to the student.")
    parser.add_argument(
        "-v",
        "--verbosity",
        dest="verbosity",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging verbosity level (choose from DEBUG, INFO, WARNING, ERROR, CRITICAL).",
    )
    parser.add_argument(
        "-w",
        "--working-dir",
        dest="working_dir",
        help="Set the working directory where the files will be created.",
    )
    parser.add_argument(
        "-a",
        "--archive-name",
        dest="archive_name",
        help="The name of the archive to be created. [default=tp4ogl-v<version>]"
    )
    parser.add_argument(
        "--skip-cleaning",
        action='store_true',
        dest="skip_cleaning",
        help="Skip the cleaning step leaving the code available.",
    )
    parser.add_argument(
        "--skip-packaging",
        action='store_true',
        dest="skip_packaging",
        help="Skip the cleaning step leaving the code available.",
    )

    args = parser.parse_args()

    if args.verbosity:
        logger.setLevel(args.verbosity)

    if args.working_dir:
        wdir = args.working_dir
    else:
        wdir = tempfile.gettempdir()

    if args.archive_name:
        zip_name = args.archive_name
    else:
        version = subprocess.check_output(["git", "describe", "--tags"]).decode("utf-8").strip()[1:]
        zip_name = f"tp4ogl-v{version}"

    try:
        main(working_dir=wdir, archive_name=zip_name, skip_cleaning=args.skip_cleaning, skip_packaging=args.skip_packaging)
    except Exception as e:
        logger.exception(e)
        sys.exit(1)
