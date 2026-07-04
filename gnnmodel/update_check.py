"""Helpers for checking app updates from GitHub Releases."""

from __future__ import annotations

import json
from dataclasses import dataclass
from urllib.request import Request, urlopen

GITHUB_REPOSITORY = "wildsonbbl/gnnepcsaftwebapp"
LATEST_RELEASE_URL = f"https://api.github.com/repos/{GITHUB_REPOSITORY}/releases/latest"


@dataclass(frozen=True)
class ReleaseInfo:
    """Minimal release metadata returned by the GitHub API."""

    tag_name: str
    html_url: str
    name: str
    body: str = ""


def _version_parts(version: str) -> tuple[int, ...]:
    cleaned = version.strip().lstrip("vV")
    parts: list[int] = []

    for part in cleaned.split("."):
        if not part.isdigit():
            break
        parts.append(int(part))

    return tuple(parts)


def is_newer_version(remote_version: str, current_version: str) -> bool:
    """Return whether the remote semantic version is newer than the current one."""

    remote_parts = _version_parts(remote_version)
    current_parts = _version_parts(current_version)

    if not remote_parts or not current_parts:
        return False

    longest_length = max(len(remote_parts), len(current_parts))
    remote_parts = remote_parts + (0,) * (longest_length - len(remote_parts))
    current_parts = current_parts + (0,) * (longest_length - len(current_parts))
    return remote_parts > current_parts


def fetch_latest_release(timeout: float = 5.0) -> ReleaseInfo:
    """Fetch the latest GitHub release for the app."""

    request = Request(
        LATEST_RELEASE_URL,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": "GNNPCSAFT-WebApp",
        },
    )

    with urlopen(request, timeout=timeout) as response:  # nosec: trusted GitHub API
        payload = json.loads(response.read().decode("utf-8"))

    return ReleaseInfo(
        tag_name=str(payload.get("tag_name", "")),
        html_url=str(payload.get("html_url", "")),
        name=str(payload.get("name", "")),
        body=str(payload.get("body", "")),
    )
